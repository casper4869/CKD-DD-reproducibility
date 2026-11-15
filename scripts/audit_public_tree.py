"""Prepublication file/credential audit. Prints paths and rule names, never values.

Run against the local release tree before upload. This is not a substitute for
human review of all public files and third-party data licenses.
"""
import hashlib
import json
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PATTERNS={
    'github_token':re.compile(r'\b(?:ghp|gho|ghu|ghs|ghr)_[A-Za-z0-9]{20,}\b'),
    'github_fine_grained_token':re.compile(r'\bgithub_pat_[A-Za-z0-9_]{30,}\b'),
    'private_key':re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
    'cloud_access_key':re.compile(r'\bAKIA[0-9A-Z]{16}\b'),
    'signed_url':re.compile(r'https?://[^\s<>\"]*[?&](?:X-Amz-Signature|X-Amz-Credential|signature|access_token|token)=[^\s<>\"]+',re.I),
    'private_absolute_path':re.compile(r'\b[A-Za-z]:[\\/]'),
    'credential_assignment':re.compile(r'''(?im)^\s*(?:password|passwd|api_key|access_token|refresh_token|client_secret)\s*[:=]\s*["'][^"']{5,}["']'''),
}
ALLOWED_SUFFIX={'.md','.txt','.csv','.json','.py','.R','.cff','.gitignore',''}


def main():
    files=[];issues=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or '.git' in p.relative_to(ROOT).parts or '__pycache__' in p.relative_to(ROOT).parts:continue
        rel=p.relative_to(ROOT).as_posix()
        if p.is_symlink():issues.append({'file':rel,'rule':'symlink'})
        if rel.startswith(('results/','.cache/')):issues.append({'file':rel,'rule':'generated_or_cached_files_in_public_tree'})
        if rel.startswith('data/private/') and rel!='data/private/README.md':issues.append({'file':rel,'rule':'private_input_in_public_tree'})
        if p.suffix not in ALLOWED_SUFFIX:issues.append({'file':rel,'rule':'unexpected_binary_or_format'})
        raw=p.read_bytes()
        try:content=raw.decode('utf-8-sig')
        except UnicodeDecodeError:
            issues.append({'file':rel,'rule':'not_utf8_text'});continue
        # Audit script contains the pattern definitions themselves; do not match
        # the definitions as actual exposed credential material.
        if rel!='scripts/audit_public_tree.py':
            for name,pattern in PATTERNS.items():
                if pattern.search(content):issues.append({'file':rel,'rule':name})
        files.append({'file':rel,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()})
    report={'status':'pass' if not issues else 'fail','file_count':len(files),'bytes':sum(x['bytes']for x in files),'issues':issues,'files':files}
    print(json.dumps(report,indent=2))
    raise SystemExit(0 if report['status']=='pass' else 1)


if __name__=='__main__':main()
