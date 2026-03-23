#!/usr/bin/env python3
import gmpy2
import shutil
import subprocess

n = 8749002899132047699790752490331099938058737706735201354674975134719667510377522805717156720453193651
e = 65537
ct = 3021569373773402689513257373362764131880473249842187164838297943840513930619586623604677697191914325

_rng = gmpy2.random_state(0xC0FFEE)

def _randrange(n):
    # return random integer in [1, n-1]
    return gmpy2.mpz_random(_rng, n - 1) + 1


def pollard_brent(n):
    if n % 2 == 0:
        return gmpy2.mpz(2)
    if n % 3 == 0:
        return gmpy2.mpz(3)

    while True:
        y = _randrange(n)
        c = _randrange(n)
        m = 128
        g = r = q = gmpy2.mpz(1)
        while g == 1:
            x = y
            for _ in range(int(r)):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(int(min(m, r - k))):
                    y = (y * y + c) % n
                    q = (q * abs(x - y)) % n
                g = gmpy2.gcd(q, n)
                k += m
            r *= 2

        if g == n:
            while True:
                ys = (ys * ys + c) % n
                g = gmpy2.gcd(abs(x - ys), n)
                if g > 1:
                    break
        if g != n:
            return g


def factor(n, out):
    if n == 1:
        return
    if gmpy2.is_prime(n):
        out.append(n)
        return
    d = pollard_brent(n)
    factor(d, out)
    factor(n // d, out)


def factor_with_gmpy2():
    factors = []
    factor(gmpy2.mpz(n), factors)
    return sorted(int(x) for x in factors)


def factor_with_sage():
    if not shutil.which("sage"):
        raise RuntimeError("sage not found")
    cmd = (
        "HOME=/tmp sage -c "
        "'n=%d; print(factor(n))'"
        % n
    )
    out = subprocess.check_output(cmd, shell=True, text=True).strip()
    # Output format: p1 * p2 * p3 * p4
    factors = [int(x.strip()) for x in out.split("*")]
    return sorted(factors)


def main():
    factors = factor_with_sage() if shutil.which("sage") else factor_with_gmpy2()

    phi = 1
    for p in factors:
        phi *= (p - 1)

    d = pow(e, -1, phi)
    plaintext = pow(ct, d, n)
    plaintext_hex = hex(plaintext)[2:]
    if len(plaintext_hex) % 2:
        plaintext_hex = "0" + plaintext_hex
    flag = bytes.fromhex(plaintext_hex).decode()

    print(f"factors: {factors}")
    print(f"flag: {flag}")


if __name__ == "__main__":
    main()
