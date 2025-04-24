from hashlib import sha3_256

# --------------------------
# ML-KEM Formulas (Core Only)
# --------------------------
n = 4  # Polynomial degree (small for testing)
q = 17  # Small modulus for illustration
k = 2  # Module rank
eta = 1  # Error range (-1, 0, 1)


def poly_add(a, b):
    return [(x + y) % q for x, y in zip(a, b)]


def poly_mul(a, b):  # a * b mod (X^n +1)
    c = [0] * n
    for i in range(n):
        for j in range(n):
            c[(i + j) % n] = (c[(i + j) % n] + a[i] * b[j]) % q
    return c


# --------------------------
# Key Generation: t = A·s + e
# --------------------------
def keygen():
    # Generate A (k x k matrix of polynomials)
    A = [
        [[1, 0, 0, 0], [2, 1, 3, 0]],  # Example static A for simplicity
        [[0, 1, 2, 0], [1, 1, 1, 1]],
    ]
    s = [[1, -1, 0, 1], [0, 1, -1, 0]]  # Secret s (k polynomials)
    e = [[1, 0, -1, 1], [0, 1, 0, -1]]  # Error e (k polynomials)

    # Compute t = A·s + e (matrix-vector multiplication)
    t = []
    for i in range(k):
        ti = [0] * n
        for j in range(k):
            aij = A[i][j]
            sj = s[j]
            ti = poly_add(ti, poly_mul(aij, sj))
        t.append(poly_add(ti, e[i]))
    return (A, t), s


# --------------------------
# Encapsulation: u = Aᵀ·r + e1, v = tᵀ·r + e2 + m_encode
# --------------------------
def encapsulate(pk):
    A, t = pk
    r = [[1, -1, 0, 0], [0, 1, -1, 0]]  # Random r
    e1 = [[0, 1, 0, -1], [1, 0, 0, 0]]  # Error e1
    m = [1, 0, 1, 0]  # Binary message

    # u = Aᵀ·r + e1 (transpose matrix multiply)
    u = []
    for i in range(k):
        ui = [0] * n
        for j in range(k):
            aji = A[j][i]  # Transpose
            rj = r[j]
            ui = poly_add(ui, poly_mul(aji, rj))
        u.append(poly_add(ui, e1[i]))

    # v = tᵀ·r + e2 + encode(m)
    v = [0] * n
    for j in range(k):
        tj = t[j]
        rj = r[j]
        v = poly_add(v, poly_mul(tj, rj))
    e2 = [0, -1, 1, 0]
    v = poly_add(v, e2)
    m_enc = [x * (q // 2) for x in m]  # Encode: 1→q/2, 0→0
    v = poly_add(v, m_enc)

    return u, v, m


# --------------------------
# Decapsulation: w = v - sᵀ·u
# --------------------------
def decapsulate(u, v, s):
    w = v.copy()
    for i in range(k):
        si = s[i]
        ui = u[i]
        w = [(wk - sij * ukj) % q for wk, sij, ukj in zip(w, si, ui)]
    # Decode: closest to q/2
    m_dec = [1 if wi > q // 4 else 0 for wi in w]
    return m_dec


# --------------------------
# Test Formula Execution
# --------------------------
pk, sk = keygen()
u, v, m = encapsulate(pk)
m_dec = decapsulate(u, v, sk)

print("Original Message:", m)
print("Decrypted Message:", m_dec)
print("Success:", m == m_dec)
