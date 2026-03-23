# ClusterRSA

## Summary

Challenge description: A message has been encrypted using RSA, but this time something feels... more crowded than usual.

This is an RSA variant where the modulus is not a product of two primes, but of four primes. Once the prime factors are recovered, the private key can be derived from Euler's totient and the ciphertext can be decrypted normally.

**Artifacts:**

- `message.txt`: provided RSA parameters
- `solve.py`: script to factor the modulus and decrypt the ciphertext

## Context

The challenge provides a standard RSA public key and ciphertext:

- `n`: the modulus
- `e`: the public exponent
- `ct`: the ciphertext

The hints imply that the modulus has more than two prime factors. This means standard RSA decryption still works, but the private key must be computed using the correct totient for a multi-prime modulus.

## Vulnerability

RSA security depends on the difficulty of factoring the modulus. In standard RSA, $n = p q$ for two large primes. Here, $n = p q r s$ with four primes. Factoring a 4-prime modulus can be easier than factoring a 2-prime modulus of similar size, because each prime is smaller and more likely to be found quickly by integer factorization algorithms.

Once the prime factors are known, the private key is computed with Euler's totient. For two primes, RSA uses:

$$
\varphi(pq) = (p-1)(q-1)
$$

This extends naturally by repeated application. Let $n_1 = pq$ and $n_2 = rs$. Then:

$$
\varphi(n_1) = (p-1)(q-1), \quad \varphi(n_2) = (r-1)(s-1)
$$

Because $n_1$ and $n_2$ are coprime, the totient is multiplicative:

$$
\varphi(n_1 n_2) = \varphi(n_1)\varphi(n_2)
$$

So for four primes:

$$
\varphi(n) = (p-1)(q-1)(r-1)(s-1)
$$

and the private exponent is:

$$
d \equiv e^{-1} \pmod{\varphi(n)}
$$

This vulnerability would fall under [CWE-327: Use of a Broken or Risky Cryptographic Algorithm](https://cwe.mitre.org/data/definitions/327.html).

## Exploitation

The following steps are implemented in [solve.py](./solve.py):

1. Factor the modulus $n$ into four primes using Pollard-Brent.

   Pollard-Brent is an optimized version of Pollard's rho factorization. It uses a pseudo-random sequence
   $x_{i+1} = x_i^2 + c \pmod{n}$ and looks for a cycle. When two values in the sequence are close, their
   difference shares a non-trivial GCD with $n$. Brent's improvement batches the GCD computations to reduce
   expensive gcd calls: it multiplies many differences together and computes a single gcd per batch.
   This is especially effective for medium-size composites like this one.

   ```python
   def pollard_brent(n):
       if n % 2 == 0:
           return 2
       if n % 3 == 0:
           return 3
       while True:
           y = _randrange(n)
           c = _randrange(n)
           m = 128
           g = r = q = 1
           while g == 1:
               x = y
               for _ in range(r):
                   y = (y * y + c) % n
               k = 0
               while k < r and g == 1:
                   ys = y
                   for _ in range(min(m, r - k)):
                       y = (y * y + c) % n
                       q = (q * abs(x - y)) % n
                   g = gcd(q, n)
                   k += m
               r *= 2
           if g == n:
               while True:
                   ys = (ys * ys + c) % n
                   g = gcd(abs(x - ys), n)
                   if g > 1:
                       break
           if g != n:
               return g
   ```
2. Compute $\varphi(n)$ as the product of $(p-1)$ for each prime factor.

   ```python
   phi = 1
   for p in factors:
       phi *= (p - 1)
   ```
3. Compute the private exponent $d = e^{-1} \bmod \varphi(n)$.

   ```python
   d = pow(e, -1, phi)
   ```
4. Decrypt the ciphertext with $m = c^d \bmod n$ and decode the resulting bytes.

   ```python
   pt = pow(ct, d, n)
   pt_hex = hex(pt)[2:]
   if len(pt_hex) % 2:
       pt_hex = "0" + pt_hex
   flag = bytes.fromhex(pt_hex).decode()
   ```

Running the script prints the factors and the recovered flag.

## Remediation

Use standard RSA with exactly two large primes or, better, use modern schemes like RSA-OAEP or elliptic-curve cryptography. Multi-prime RSA can be safe with proper key sizing, but for CTF-sized parameters it is significantly easier to factor.
