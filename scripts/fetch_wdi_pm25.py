"""Optional, opt-in World Bank API access; no crawlers or GBD downloading.

Only official API requests are supported. Existing cache works offline.
No component of the analysis pipeline calls this module.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode, urlparse
from urllib.request import Request, build_opener, HTTPRedirectHandler

API = 'https://api.worldbank.org/v2/country/all/indicator/EN.ATM.PM25.MC.M3'
MIN_INTERVAL = 10.0


def require_official_url(url):
    u = urlparse(url)
    if u.scheme != 'https' or u.hostname != 'api.worldbank.org' or u.username or u.password:
        raise ValueError('Only the official HTTPS World Bank API is permitted.')


class OfficialRedirectsOnly(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        require_official_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def retry_after_seconds(value, fallback):
    if not value:
        return fallback
    try:
        return max(0.0, float(value))
    except ValueError:
        try:
            return max(0.0, (parsedate_to_datetime(value) - datetime.now(timezone.utc)).total_seconds())
        except (ValueError, TypeError, OverflowError):
            return fallback


class Client:
    def __init__(self, cache, online=False, interval=MIN_INTERVAL, opener=None, sleep=time.sleep, clock=time.monotonic):
        if interval < MIN_INTERVAL:
            raise ValueError('Requests must be at least 10 seconds apart.')
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.online, self.interval, self.last = online, interval, None
        self.opener = opener or build_opener(OfficialRedirectsOnly())
        self.sleep, self.clock = sleep, clock

    def get(self, url):
        require_official_url(url)
        path = self.cache / (hashlib.sha256(url.encode()).hexdigest() + '.json')
        if path.exists():
            return json.loads(path.read_text(encoding='utf-8'))['response']
        if not self.online:
            raise RuntimeError('Offline mode: response is not cached. Explicit --online is required.')
        for attempt in range(4):
            if self.last is not None:
                self.sleep(max(0.0, self.interval - (self.clock() - self.last)))
            self.last = self.clock()
            try:
                request = Request(url, headers={'Accept': 'application/json', 'User-Agent': 'CKD-DD-reproducibility/1.0 (official-API-only; no crawler)'})
                with self.opener.open(request, timeout=90) as response:
                    require_official_url(response.geturl())
                    data = json.loads(response.read())
                tmp = path.with_suffix('.tmp')
                tmp.write_text(json.dumps({'url': url, 'retrieved_utc': datetime.now(timezone.utc).isoformat(), 'response': data}), encoding='utf-8')
                tmp.replace(path)
                return data
            except HTTPError as exc:
                if exc.code not in (429, 500, 502, 503, 504) or attempt == 3:
                    raise
                delay = max(self.interval, retry_after_seconds(exc.headers.get('Retry-After'), 60.0 * 2**attempt))
                self.sleep(delay)
        raise RuntimeError('Retry budget exhausted.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--online', action='store_true', help='Explicitly allow serial official-API requests. Follow all service terms and rate limits. No crawlers.')
    p.add_argument('--cache', type=Path, default=Path('.cache/wdi'))
    p.add_argument('--output', type=Path, default=Path('results/wdi_current_pm25.json'))
    args = p.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    lock = args.cache / '.serial-request.lock'
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise SystemExit('Another request process or stale lock exists. Do not run requests concurrently.')
    try:
        os.close(fd)
        client = Client(args.cache, online=args.online)
        all_rows, pages = [], 1
        page = 1
        while page <= pages:
            url = API + '?' + urlencode({'date': '1992:2020', 'format': 'json', 'per_page': '20000', 'page': page})
            data = client.get(url)
            if not isinstance(data, list) or len(data) != 2 or not isinstance(data[1], list):
                raise RuntimeError('Unexpected API response; stop without further requests.')
            pages = int(data[0]['pages'])
            if not 1 <= pages <= 10:
                raise RuntimeError('Unexpected page count; stop without crawling.')
            all_rows.extend(data[1])
            page += 1
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps({'warning': 'Current WDI API data, not the frozen study PM2.5 snapshot. No study input is overwritten.', 'indicator': 'EN.ATM.PM25.MC.M3', 'unit': 'micrograms per cubic metre', 'records': all_rows}, indent=2), encoding='utf-8')
    finally:
        lock.unlink(missing_ok=True)


if __name__ == '__main__':
    main()
