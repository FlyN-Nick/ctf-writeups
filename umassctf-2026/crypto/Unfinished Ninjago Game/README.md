# Unfinished Ninjago Game

## Summary

I'm writing a ninjago game, but apparently there's some crazy cheat codes someone has already found!? And apparently that's supposed to help decrypt hidden secret stuff by the Overlord?

**Artifacts:**

- `main`: ELF 64-bit PIE x86-64 binary; the challenge server
- `PRNG.h`: source for the modified xoshiro512\*\* PRNG embedded in `main`
- `exploit.py`: Python script that implements the attack

## Context

### PRNG

`PRNG.h` implements a modified [xoshiro512\*\*](https://prng.di.unimi.it/) generator. The state is eight 64-bit words `s[0..7]` (a total of 512 bits), seeded by `getrandom` before anything else happens. The file also defines thirty-two precomputed **jump functions** (`jump0`, `jump16`, …, `jump496`) that advance the state by fixed powers of two without calling `next()` iteratively. The comment in the file hints at a planned-but-unimplemented feature:

```c
// TODO: write short code to compute pow(x, n, char_poly(next)) for arbitrary n
```

This is a red herring; the intended exploit does not require computing arbitrary jumps.

### Program Flow

On each TCP connection, `main` runs the following sequence:

1. **Seed** — `getrandom(s, 64, 0)` fills all eight state words with fresh random bytes.
2. **`load()`** — Opens `flag.txt`, reads up to 64 bytes into a local buffer, and prints the ciphertext:

   ```c
   for (int i = 0; i < 8; i++)
       printf("%llu ", local_buf[i] ^ s[i]);
   putc('\n', stdout);
   remove("flag.txt");
   ```

   The flag is treated as eight little-endian `uint64` words. The PRNG state is **never advanced** before or after this step, so `s[i]` at this point is exactly the seeded value.
3. **`explore_right()`** — An interactive command loop that reads single characters from stdin:

   | Input         | Action                                                                                     |
   | ------------- | ------------------------------------------------------------------------------------------ |
   | `e`         | exit                                                                                       |
   | `m`         | call `explore_middle()` — the oracle                                                    |
   | `a`         | advance PRNG with a**randomly chosen** jump (via `getrandom`; not user-controlled) |
   | `s`         | set a global pointer `d` to a stack address                                              |
   | `w`         | write one `getc` byte via `d` (out-of-bounds write, but cannot reach the canary)       |
   | `r` / `l` | recurse into `explore_right` / `explore_left`                                          |

### The Oracle

`explore_middle()` computes and outputs a single byte:

```c
uint64_t sum = s[0] + s[1] + ... + s[7];   // wrapping uint64 addition
putc((int)(sum % 101), stdout);
```

The magic-number division constant `0x446f86562d9faee5` is the compiler's fast unsigned divide-by-101. If `'m'` is sent **before** any `'a'` commands, the oracle reflects the original seeded state — the same state used to encrypt the flag.

## Vulnerability

Because the attacker knows `ct[i] = s[i] XOR flag_word[i]` and the oracle $O = \left(\sum_{i=0}^{7} s_i\right) \bmod 101$, all the information needed to recover the flag is present in each connection — it just needs to be extracted algebraically.

### Bit-level decomposition

Index the 512 flag bits as $j = 0, 1, \ldots, 511$, where bit $j$ belongs to word $\lfloor j / 64 \rfloor$ at bit position $j \bmod 64$. Let:

- $b_j = \text{bit}_j(ct)$ — the corresponding ciphertext bit (known)
- $x_j = \text{bit}_j(flag)$ — the flag bit (unknown)

Since XOR is bitwise addition mod 2, $\text{bit}_j(s) = b_j \oplus x_j$. As integers, this is:

$$
b_j \oplus x_j = b_j + x_j - 2 b_j x_j
$$

Substituting into the oracle (where bit $j$ contributes weight $2^{j \bmod 64}$ to its containing word):

$$
O = \sum_{j=0}^{511} 2^{j \bmod 64} \cdot (b_j + x_j - 2 b_j x_j) \pmod{101}
$$

$$
O - \underbrace{\sum_{j=0}^{511} 2^{j \bmod 64} \cdot b_j}_{\displaystyle=\;\sum_{i=0}^{7} ct_i \pmod{101}} \;\equiv\; \sum_{j=0}^{511} \underbrace{(1 - 2b_j) \cdot 2^{j \bmod 64}}_{\displaystyle H_j} \cdot x_j \pmod{101}
$$

### Linear system

Defining:

$$
h \;=\; \left(O - \sum_{i=0}^{7} ct_i\right) \bmod 101 \qquad H_j \;=\; (1 - 2b_j) \cdot 2^{j \bmod 64} \bmod 101
$$

every connection yields one linear equation $\sum_j H_j x_j \equiv h \pmod{101}$, with known coefficients determined by that connection's ciphertext and oracle reading. Because $101$ is prime, $\mathbb{Z}/101\mathbb{Z}$ is a field. With 512 fresh connections producing linearly independent equations (guaranteed with overwhelming probability because the ciphertexts are independently random), the $512 \times 512$ system has a **unique solution** over $\mathbb{Z}/101\mathbb{Z}$. Since the unknowns $x_j$ are flag bits, they are in $\{0, 1\} \subset \mathbb{Z}/101\mathbb{Z}$, so the unique field solution equals the true bit values.

## Exploitation

The below steps are implemented in [exploit.py](./exploit.py).

### 1. Collect (ciphertext, oracle) pairs

Open 20 parallel TCP connections. For each connection, read the eight ciphertext words and send `'m'` immediately (before any `'a'` commands so the state is still at its initial seeded value), then read the oracle byte.

```python
def query():
    """Open one connection, return (ct, oracle) or raise on error."""
    s = socket.socket()
    s.settimeout(15)
    try:
        s.connect((HOST, PORT))
        data = b""
        while b"\n" not in data:
            data += s.recv(4096)
        ct = [int(x) for x in data.split(b"\n")[0].split()]
        s.sendall(b"m")
        oracle = s.recv(1)[0]
        s.sendall(b"e")
        return ct, oracle
    finally:
        s.close()
```

Worker threads run `query()` in a loop and push results onto a shared queue. The main thread reads from the queue until 560 samples are collected (48 extra beyond 512 as a safety margin).

### 2. Build the linear system

Translate each `(ct, oracle)` pair into a row of the matrix $H$ and the right-hand-side vector $h$ over $\mathbb{Z}/101\mathbb{Z}$.

The ciphertext values can exceed $2^{63}$, so Python's arbitrary-precision integers are used for the modular arithmetic in $h$ rather than NumPy (which would silently overflow `int64`).

```python
def build_system(samples):
    """
    Build H (n×512) and h (n,) over Z/101Z from (ct, oracle) pairs.

    For flag bit j spanning word index j//64 at bit position j%64:
      H[i][j] = (1 - 2 * bit_j(ct[i])) * 2^(j%64)  mod 101
      h[i]    = (oracle[i] - Σ_k ct[i][k])           mod 101
    """
    n = len(samples)
    H = np.zeros((n, N_BITS), dtype=np.int64)
    h = np.zeros(n, dtype=np.int64)
    for i, (ct, oracle) in enumerate(samples):
        h[i] = (oracle - sum(c % MOD for c in ct)) % MOD  # Python int: safe for uint64
        for j in range(N_BITS):
            ct_bit = (ct[j // 64] >> (j % 64)) & 1
            H[i, j] = (1 - 2 * ct_bit) * POW2[j] % MOD
    return H, h
```

`POW2[j] = pow(2, j % 64, 101)` is precomputed once at module load.

### 3. Solve via RREF over $\mathbb{Z}/101\mathbb{Z}$

Gaussian elimination (reduced row echelon form) is performed on the augmented matrix $[H \mid h]$. For each pivot column, the pivot row is scaled so the leading entry becomes 1 (using Fermat's little theorem to compute the modular inverse: $a^{-1} \equiv a^{99} \pmod{101}$), and then that column is eliminated from every other row simultaneously using NumPy vectorised operations.

```python
def rref_mod_p(A, b, p):
    """
    Reduce augmented matrix [A|b] to RREF over Z/pZ (p prime).
    Returns (reduced matrix, list of pivot column indices).
    """
    M = np.hstack([A % p, (b % p).reshape(-1, 1)]).astype(np.int64)
    pivots, r = [], 0

    for col in range(A.shape[1]):
        # Find a non-zero entry in this column at or below row r
        nz = np.nonzero(M[r:, col] % p)[0]
        if not len(nz):
            continue

        M[[r, r + nz[0]]] = M[[r + nz[0], r]]                 # swap into pivot position
        M[r] = M[r] * pow(int(M[r, col]), p - 2, p) % p       # scale so pivot = 1

        # Eliminate this column in every other row
        mask = (M[:, col] % p != 0)
        mask[r] = False
        M[mask] = (M[mask] - M[mask, col:col+1] * M[r]) % p

        pivots.append(col)
        r += 1
        if r >= M.shape[0]:
            break

    return M, pivots
```

After RREF the solution vector is read directly from the last column of the reduced matrix: `x[col] = M[i, 512]` for each pivot `(i, col)`.

### 4. Reconstruct the flag

The 512 recovered bits are packed back into 64 bytes. Bit $j$ lands at byte $\lfloor j/8 \rfloor$ at bit position $j \bmod 8$, which follows directly from the little-endian layout of the eight `uint64` flag words.

```python
def reconstruct_flag(x):
    """
    Pack the 512 recovered flag bits back into 64 bytes.
    Bit j maps to byte j//8 (little-endian within each 64-bit word), bit position j%8.
    """
    flag = bytearray(64)
    for j in range(N_BITS):
        flag[j // 8] |= x[j] << (j % 8)
    return bytes(flag)
```

Running the exploit against the live server:

```
[*] Collecting 560 samples (20 threads)...
    100/560  (221 q/s)
    200/560  (225 q/s)
    ...
[*] Collected 560 samples in 2.5s
[*] Building 560×512 matrix over Z/101Z...
[*] Solving via RREF...

[+] UMASS{sparse_fourier_transforms_are_so_much_fun!fhwtftw!yayayay}
```

## Remediation

The root cause is that the oracle `(Σ s[i]) mod 101` leaks a linear combination of the PRNG state, and because each connection re-seeds the PRNG independently while broadcasting the ciphertext, an attacker accumulates enough linear equations to fully determine the secret state. Two concrete fixes:

1. **Use a larger modulus.** A modulus of 101 means coefficients lie in a tiny field, making the linear system easy to solve. If the modulus were a 64-bit prime (or the output were a full 64-bit word), the attacker would need far more information per query and the system would be much harder to solve in practice.
2. **Do not expose a linear oracle over freshly seeded states with known ciphertext.** The combination of a known linear function of the state *and* a known linear function of `state XOR flag` is what enables the attack. Either hide one of the two (don't print the ciphertext, or don't offer the oracle), or introduce non-linearity (e.g., hash the state before outputting it).
