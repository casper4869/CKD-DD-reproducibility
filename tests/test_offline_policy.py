"""Deterministic tests of optional API safeguards; no network requests."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from urllib.error import HTTPError

spec=importlib.util.spec_from_file_location('wdi',Path(__file__).resolve().parents[1]/'scripts/fetch_wdi_pm25.py')
wdi=importlib.util.module_from_spec(spec);spec.loader.exec_module(wdi)


class Response:
    def __enter__(self):return self
    def __exit__(self,*args):pass
    def read(self):return json.dumps([{'pages':1},[]]).encode()
    def geturl(self):return wdi.API


class Opener:
    def __init__(self,rate_limit=False):self.calls=0;self.rate_limit=rate_limit
    def open(self,req,timeout):
        self.calls+=1
        if self.rate_limit and self.calls==1:
            raise HTTPError(req.full_url,429,'limited',{'Retry-After':'45'},None)
        return Response()


class PolicyTests(unittest.TestCase):
    def test_default_offline_makes_no_requests(self):
        with tempfile.TemporaryDirectory() as d:
            op=Opener();client=wdi.Client(d,opener=op)
            with self.assertRaises(RuntimeError):client.get(wdi.API)
            self.assertEqual(op.calls,0)
    def test_official_https_domain_only(self):
        for url in ['http://api.worldbank.org/v2','https://evil.example/v2','https://api.worldbank.org.evil.example/v2','https://user:password@api.worldbank.org/v2']:
            with self.assertRaises(ValueError):wdi.require_official_url(url)
    def test_cache_and_pacing(self):
        with tempfile.TemporaryDirectory() as d:
            op=Opener();waits=[];client=wdi.Client(d,online=True,opener=op,sleep=waits.append,clock=lambda:0)
            client.get(wdi.API+'?page=1');client.get(wdi.API+'?page=1')
            self.assertEqual(op.calls,1)
            client.get(wdi.API+'?page=2');self.assertEqual(op.calls,2)
            self.assertEqual(waits,[10.0])
            self.assertEqual(wdi.Client(d,opener=Opener()).get(wdi.API+'?page=1'),[{'pages':1},[]])
    def test_rate_limit_honors_retry_after(self):
        with tempfile.TemporaryDirectory() as d:
            waits=[];op=Opener(True);client=wdi.Client(d,online=True,opener=op,sleep=waits.append,clock=lambda:0)
            client.get(wdi.API)
            self.assertEqual(op.calls,2);self.assertGreaterEqual(max(waits),45)
    def test_cannot_reduce_minimum_interval(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):wdi.Client(d,interval=0)


if __name__=='__main__':unittest.main()
