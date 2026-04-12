# Hens and Roosters

## Summary

The challenge presents a web application where users must accumulate seven "studs" to purchase a Lego set and obtain the flag. Studs are earned by submitting valid signatures for the current stud count. While the server provides free signatures for the first three studs (0, 1, and 2), reaching seven requires bypassing rate limits, exploiting a race condition, and leveraging a cryptographic property of the signature scheme to "clone" valid signatures.

**Artifacts:**

- `backend/app.py`: Flask application handling user progress and signature verification.
- `backend/uov.py`: Implementation of the Unbalanced Oil and Vinegar (UOV) signature scheme.
- `proxy/haproxy.cfg`: HAProxy configuration with a vulnerable rate-limiting rule.
- `solve.py`: Exploit script demonstrating the full attack chain.

## Context

The application uses the **Unbalanced Oil and Vinegar (UOV)** signature scheme over the extension field $\mathbb{F}_{2^7}$. A user is identified by a `uid`, and their progress (number of studs) is stored in Redis. 

To gain a stud, a user must `POST` a valid signature for the payload `str(studs) + '|' + uid` to the `/work` endpoint. The server implements a caching mechanism in Redis to speed up verification of recently seen signatures:

```python
@app.post('/work')
def work():
    # ... read studs and payload ...
    value = r.get(str(sig))
    if value is None:
        r.set(sig, b'-', ex=240)
        verified = uov.verify(payload, sig_bytes) # Slow path (~2.5s)
        if verified:
            r.set(sig, payload, ex=240)
    else:
        verified = value.decode() == payload # Fast path (cache hit)
    
    if verified:
        studs = r.incr(uid)
        # ... return next free signature if studs <= 2 ...
```

The server is protected by HAProxy, which enforces a rate limit of one request every 20 seconds.

## Vulnerability

Three distinct vulnerabilities are chained to achieve the exploit.

### 1. HAProxy Rate-Limit Bypass - [CWE-837: Improper Enforcement of a Single, Unique Action](https://cwe.mitre.org/data/definitions/837.html)

The HAProxy configuration tracks request rates based on the full URL, including query parameters:

```haproxy
stick-table type string len 2048 size 100k expire 20s store http_req_rate(20s)
http-request track-sc0 url
http-request deny deny_status 429 if { sc_http_req_rate(0) gt 1 }
```

By appending unique query parameters (e.g., `/work?x=1`, `/work?x=2`), an attacker can force HAProxy to treat each request as a distinct URL, effectively bypassing the rate limit and enabling high-concurrency attacks.

### 2. TOCTOU Race Condition in `/work` - [CWE-367: Time-of-Check Time-of-Use (TOCTOU) Race Condition](https://cwe.mitre.org/data/definitions/367.html)

The `/work` endpoint reads the current stud count from Redis at the beginning of the request and only increments it after a successful (and slow) signature verification. 

```python
studs = r.get(uid)                      # (1) Read current count
payload = str(studs) + '|' + uid
# ... slow uov.verify(payload, sig) ... # (2) Large window (~2.5s)
if verified:
    studs = r.incr(uid)                 # (3) Increment count
```

Because verification takes several seconds, multiple concurrent requests can all read the same `studs` value (e.g., `2`), verify their respective signatures for the same payload (`"2|uid"`), and then all trigger the increment. This allows jumping from 2 studs to 7 studs in a single race window.

### 3. UOV Frobenius Signature Cloning - [CWE-327: Use of a Broken or Risky Cryptographic Algorithm](https://cwe.mitre.org/data/definitions/327.html)

The UOV implementation uses public key matrices $P_i$ with coefficients restricted to the base field $\mathbb{F}_2 \subseteq \mathbb{F}_{2^7}$. Additionally, the verification targets $t_i$ are bits derived from a hash, so $t_i \in \{0, 1\}$.

In any field of characteristic 2, the Frobenius automorphism $\sigma(x) = x^2$ is linear. If a signature $\mathbf{x}$ satisfies the verification equation $\mathbf{x}^\top P_i \mathbf{x} = t_i$, then applying $\sigma$ to $\mathbf{x}$ yields:

$$\sigma(\mathbf{x})^\top P_i \sigma(\mathbf{x}) = \sigma(\mathbf{x}^\top P_i \mathbf{x}) = \sigma(t_i) = t_i$$

The equality $\sigma(P_i) = P_i$ and $\sigma(t_i) = t_i$ holds because their elements are in $\mathbb{F}_2$, which is fixed by $\sigma$. Since $[\mathbb{F}_{2^7} : \mathbb{F}_2] = 7$, an attacker can generate 6 additional valid signatures for the *same* message from a single observed signature. These "cloned" signatures are not in the server's Redis cache, forcing every concurrent request into the slow verification path and widening the race window.

## Exploitation

The exploit is implemented in `solve.py` and follows these steps:

1.  **Preparation**: Use the rate-limit bypass to quickly obtain `sig_2` (the signature for `"2|uid"`) from the server by advancing from 0 to 2 studs.
2.  **Signature Cloning**: Compute the Frobenius orbit of `sig_2` by repeatedly squaring each byte in $\mathbb{F}_{2^7}$ using the irreducible polynomial $x^7 + x + 1$. This produces 5 new, uncached signatures for the payload `"2|uid"`.
3.  **The Race**: Submit the 5 cloned signatures concurrently using unique query parameters. 
    - Each request bypasses HAProxy's rate limit.
    - Each request hits the "slow path" in `/work` because the cloned signatures are uncached.
    - All 5 requests read `studs = 2` from Redis before any of them completes verification.
4.  **Completion**: Once the verifications finish, all 5 requests call `r.incr(uid)`, advancing the stud count to 7. The attacker then calls `/buy` to retrieve the flag.

```python
def _gf128_sq(a: int) -> int:
    """Square an element of GF(2^7) under x^7 + x + 1."""
    p, b = 0, a
    for _ in range(7):
        if b & 1: p ^= a
        b >>= 1
        hi = a >> 6
        a = (a << 1) & 0x7F
        if hi: a ^= 0x03
    return p

def frob_clone(sig_hex: str):
    """Generate a new valid signature using Frobenius."""
    return bytes(_gf128_sq(b) for b in bytes.fromhex(sig_hex)).hex()
```

## Remediation

1.  **Secure Rate Limiting**: Configure HAProxy to track request rates by `path` or `src` (IP address) rather than the full `url` to prevent query-parameter-based bypasses.
2.  **Atomic State Updates**: Use atomic Redis operations or Lua scripts to ensure that the "read-verify-increment" cycle is protected against race conditions. For example, use a lock or check that the stud count hasn't changed before incrementing.
3.  **Proper UOV Parameterization**: Ensure that the secret linear transformation $T$ (and consequently the public key) is chosen with entries from the full extension field $\mathbb{F}_{2^7}$ rather than being restricted to the subfield $\mathbb{F}_2$. This breaks the Frobenius symmetry.
