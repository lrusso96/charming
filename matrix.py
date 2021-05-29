from charm.toolbox.pairinggroup import ZR, G1, G2, GT, pair

from collections import namedtuple
from functools import reduce


MM_GROUP = namedtuple('MM_GROUP', ['G', 'g1', 'g2', 'gT'])


class MatrixDistribution():

    def sample(self, gtype=ZR) -> list:
        raise NotImplementedError


class DK_MDDH(MatrixDistribution):

    def __init__(self, k: int, group: MM_GROUP):
        self.k = k
        self.group = group
        self.gtype = ZR

    def sample(self) -> list:
        return [[self.group.G.random(self.gtype) for i in range(self.k)] for j in range(self.k+1)]


class MM():

    def __init__(self, group: MM_GROUP):
        self.group = group

    def __sample(self, n, gtype=ZR):
        return [self.group.G.random(gtype) for i in range(n)]

    def sample(self, rows: int, cols: int = 1, gtype=ZR):
        return MMMatrix(self.group, [self.__sample(cols, gtype) for i in range(rows)], gtype)

    def sample_from(self, distribution: MatrixDistribution):
        return MMMatrix(self.group, distribution.sample(), distribution.gtype)

    def T(M):
        return [[row[j] for row in M] for j in range(len(M[0]))]

    def get_generator(self, gtype):
        if gtype == G1:
            return self.group.g1
        elif gtype == G2:
            return self.group.g2
        elif gtype == GT:
            return self.group.gT

    def exponentiate_vector(v, base):
        return list(map(lambda x: base ** x, v))

    def exponentiate_matrix(M, base):
        return [MM.exponentiate_vector(row, base) for row in M]

    def pair_vector(v, base):
        return list(map(lambda x: pair(x, base), v))

    def pair_matrix(M, base):
        return [MM.pair_vector(row, base) for row in M]

    def rowvec_mul_mat(v, M) -> list:
        return [sum([v[i]*M[i][j] for i in range(len(v))]) for j in range(len(M[0]))]

    def rowvec_pair_mat(v, M) -> list:
        return [reduce(lambda a, b: a*b, [pair(v[i], M[i][j]) for i in range(len(v))]) for j in range(len(M[0]))]

    def rowvec_exp_mat(v, M) -> list:
        return [reduce(lambda a, b: a*b, [v[i]**M[i][j] for i in range(len(v))]) for j in range(len(M[0]))]

    def mat_mul_mat(A, B):
        return [MM.rowvec_mul_mat(row, B) for row in A]

    def mat_exp_mat(A, B):
        return [MM.rowvec_exp_mat(row, B) for row in A]

    def mat_pair_mat(A, B):
        return [MM.rowvec_pair_mat(row, B) for row in A]

    def mat_add_mat(A, B):
        return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

    def mat_sub_mat(A, B):
        return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

    def mat_had_mat(A, B):
        return [[A[i][j] * B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

    def mat_had_div_mat(A, B):
        return [[A[i][j] / B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


class MMMatrix(MM):
    def __init__(self, group: MM_GROUP, M, gtype):
        super().__init__(group)
        self.M = M
        self.gtype = gtype

    def T(self):
        return MMMatrix(self.group, MM.T(self.M), self.gtype)

    def __or__(self, o):
        if len(o.M[0]) == len(self.M[0]):
            return MMMatrix(self.group, self.M + o.M, self.gtype)
        raise NotImplementedError

    def __rshift__(self, gtype):
        if self.gtype == ZR and gtype in (G1, G2, GT):
            return MMMatrix(self.group, MM.exponentiate_matrix(self.M, self.get_generator(gtype)), gtype)
        if self.gtype == G1 and gtype == GT:
            return MMMatrix(self.group, MM.pair_matrix(self.M, self.group.g2), GT)
        if self.gtype == G2 and gtype == GT:
            return MMMatrix(self.group, MM.pair_matrix(self.M, self.group.g1), GT)
        raise NotImplementedError

    def __eq__(self, o) -> bool:
        if type(o) == MMMatrix:
            return self.gtype == o.gtype and self.M == o.M
        return False

    def __add__(self, o):
        if type(o) == MMMatrix:
            if self.gtype == o.gtype:
                if o.gtype == ZR:
                    return MMMatrix(self.group, MM.mat_add_mat(self.M, o.M), ZR)
                return MMMatrix(self.group, MM.mat_had_mat(self.M, o.M), self.gtype)
        raise NotImplementedError

    def __sub__(self, o):
        if type(o) == MMMatrix:
            if self.gtype == o.gtype:
                if o.gtype == ZR:
                    return MMMatrix(self.group, MM.mat_sub_mat(self.M, o.M), ZR)
                return MMMatrix(self.group, MM.mat_had_div_mat(self.M, o.M), self.gtype)
        raise NotImplementedError

    def __mul__(self, o):
        if type(o) == MMMatrix:
            if self.gtype == o.gtype == ZR:
                return MMMatrix(self.group, MM.mat_mul_mat(self.M, o.M), ZR)
            elif self.gtype != ZR and o.gtype == ZR:
                return MMMatrix(self.group, MM.mat_exp_mat(self.M, o.M), self.gtype)
            elif self.gtype == ZR and o.gtype != ZR:
                return (o.T() * self.T()).T()
            elif (self.gtype, o.gtype) in ((G1, G2), (G2, G1)):
                return MMMatrix(self.group, MM.mat_pair_mat(self.M, o.M), GT)
        raise NotImplementedError

    def pair_with(self, o):
        if type(o) == MMMatrix:
            if set((self.gtype, o.gtype)) == set((G1, G2)):
                return self.T() * o
        raise NotImplementedError

    def __len__(self) -> int:
        return len(self.M)

    def __str__(self) -> str:
        return str(self.M)
