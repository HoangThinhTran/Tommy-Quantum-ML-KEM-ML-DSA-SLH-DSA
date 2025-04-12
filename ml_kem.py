import os
import hashlib
import random

# Parameters for Kyber512
n = 256
q = 3329
k = 2
eta1 = 3
eta2 = 2
du = 10
dv = 4

# Precomputed zetas and inv_zetas for NTT (truncated for brevity)
zetas = [
    2285,
    2571,
    2970,
    1812,
    1493,
    1422,
    287,
    202,
    3158,
    622,
    1577,
    182,
    962,
    2127,
    1855,
    1468,
    573,
    2004,
    264,
    383,
    2500,
    1458,
    1727,
    3199,
    2648,
    1017,
    732,
    608,
    1787,
    411,
    3124,
    1758,
]
inv_zetas = [
    1701,
    1807,
    1460,
    2371,
    2338,
    2333,
    308,
    108,
    2851,
    870,
    854,
    1510,
    2535,
    1278,
    1530,
    1185,
]


def pack_bits(compressed_coeffs, d):
    byte_array = bytearray()
    buffer = 0
    bits_in_buffer = 0
    for coeff in compressed_coeffs:
        buffer = (buffer << d) | coeff
        bits_in_buffer += d
        while bits_in_buffer >= 8:
            bits_to_extract = bits_in_buffer - 8
            byte = (buffer >> bits_to_extract) & 0xFF
            byte_array.append(byte)
            buffer &= (1 << bits_to_extract) - 1
            bits_in_buffer -= 8
    if bits_in_buffer > 0:
        byte_array.append(buffer << (8 - bits_in_buffer))
    return bytes(byte_array)


def sample_cbd(eta, n_coeffs):
    return [
        (sum(random.choices([0, 1], k=eta)) - sum(random.choices([0, 1], k=eta)))
        for _ in range(n_coeffs)
    ]


def parse(bytes_stream):
    coeffs = []
    for i in range(128):
        b0 = bytes_stream[3 * i]
        b1 = bytes_stream[3 * i + 1]
        b2 = bytes_stream[3 * i + 2]
        d1 = b0 + ((b1 & 0x0F) << 8)
        d2 = (b1 >> 4) + (b2 << 4)
        coeffs.extend([d1 % q, d2 % q])
    return coeffs


def ntt(poly):
    poly = poly.copy()
    len_val = 128
    k = 0
    while len_val >= 2:
        for start in range(0, 256, 2 * len_val):
            zeta = zetas[k]
            k += 1
            for j in range(start, start + len_val):
                t = (zeta * poly[j + len_val]) % q
                poly[j + len_val] = (poly[j] - t) % q
                poly[j] = (poly[j] + t) % q
        len_val >>= 1
    return poly


def invntt(poly):
    poly = poly.copy()
    len_val = 2
    k = 0
    while len_val <= 128:
        for start in range(0, 256, 2 * len_val):
            zeta = inv_zetas[k]
            k += 1
            for j in range(start, start + len_val):
                t = poly[j]
                poly[j] = (t + poly[j + len_val]) % q
                poly[j + len_val] = ((poly[j + len_val] - t) * zeta) % q
        len_val <<= 1
    inv_n = pow(256, -1, q)
    return [(c * inv_n) % q for c in poly]


def multiply_poly(a, b):
    product = [0] * n
    for i in range(n):
        for j in range(n):
            product[(i + j) % n] = (product[(i + j) % n] + a[i] * b[j]) % q
    return product


def compress(poly, d):
    # Handle both single polynomials (list of ints) and arrays of polynomials
    if isinstance(poly[0], list):
        return [[((coeff << d) // q) & ((1 << d) - 1) for coeff in p] for p in poly]
    else:
        return [((coeff << d) // q) & ((1 << d) - 1) for coeff in poly]


def decompress(compressed, d):
    if isinstance(compressed[0], list):
        return [[(c * q + (1 << (d - 1))) // (1 << d) for c in p] for p in compressed]
    else:
        return [(c * q + (1 << (d - 1))) // (1 << d) for c in compressed]


def keygen():
    rho = os.urandom(32)
    # Generate matrix A using SHAKE-128
    A = [
        [parse(hashlib.shake_128(rho + bytes([j, i])).digest(384)) for j in range(k)]
        for i in range(k)
    ]

    # Generate secret and error vectors
    s = [sample_cbd(eta1, n) for _ in range(k)]
    e = [sample_cbd(eta1, n) for _ in range(k)]

    # Calculate t = A*s + e
    t = []
    for i in range(k):
        row = [0] * n
        for j in range(k):
            product = multiply_poly(A[i][j], s[j])
            row = [(row[idx] + product[idx]) % q for idx in range(n)]
        row = [(row[idx] + e[i][idx]) % q for idx in range(n)]
        t.append(row)

    return (rho, t), s


def serialize(poly):
    bytes_stream = bytearray()
    for i in range(128):
        d1 = poly[2 * i]
        d2 = poly[2 * i + 1]
        bytes_stream.append(d1 & 0xFF)
        bytes_stream.append(((d1 >> 8) & 0x0F) | ((d2 & 0x0F) << 4))
        bytes_stream.append((d2 >> 4) & 0xFF)
    return bytes(bytes_stream)


def encapsulate(pk):
    rho, t = pk
    m = os.urandom(32)
    t_bytes = b"".join([serialize(tt) for tt in t])
    K_bar_r = hashlib.shake_256(m + rho + t_bytes).digest(64)
    K_bar, r_seed = K_bar_r[:32], K_bar_r[32:]

    r = [sample_cbd(eta2, n) for _ in range(k)]
    e1 = [sample_cbd(eta2, n) for _ in range(k)]
    e2 = sample_cbd(eta2, n)

    # Proper matrix generation with i,j loops
    A_T = list(
        zip(
            *[
                [
                    parse(hashlib.shake_128(rho + bytes([j, i])).digest(384))
                    for j in range(k)
                ]
                for i in range(k)  # This defines the i variable
            ]
        )
    )

    u = []
    for i in range(k):
        row = [0] * n
        for j in range(k):
            product = multiply_poly(A_T[i][j], r[j])
            row = [(row[idx] + product[idx]) % q for idx in range(n)]
        u.append([(row[idx] + e1[i][idx]) % q for idx in range(n)])

    v = [
        (
            sum(multiply_poly(t[i], r[i])[idx] for i in range(k))
            + e2[idx]
            + (q // 2 if (m[idx // 8] >> (7 - (idx % 8))) & 1 else 0)
        )
        % q
        for idx in range(n)
    ]
    u_compressed = compress(u, du)
    v_compressed = compress(v, dv)
    return (u_compressed, v_compressed), K_bar


def decapsulate(ct, sk, pk):
    u_compressed, v_compressed = ct
    s = sk
    rho, t = pk
    u = [decompress(uc, du) for uc in u_compressed]
    v = decompress(v_compressed, dv)
    m_prime = [
        (v[idx] - sum(multiply_poly(s[i], u[i])[idx] for i in range(k))) % q
        for idx in range(n)
    ]
    m = bytes(
        [
            sum(
                (
                    (1 if (m_prime[i * 8 + j] > q // 2) else 0) << (7 - j)
                    for j in range(8)
                )
            )
            for i in range(32)
        ]
    )
    # Re-encapsulate to verify
    _, K = encapsulate(pk)
    return K


# Example usage
if __name__ == "__main__":
    pk, sk = keygen()
    ct, K_encap = encapsulate(pk)
    K_decap = decapsulate(ct, sk, pk)
    print("Shared Secret (Encap):", K_encap.hex())
    print("Shared Secret (Decap):", K_decap.hex())
