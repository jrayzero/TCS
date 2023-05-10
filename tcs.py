class Params(object):

    def __init__(self, N, perms, origin, refinement, coarsening, extents):
        self.N = N
        self.perms = perms
        self.origin = origin
        self.refinement = refinement
        self.coarsening = coarsening
        self.extents = extents

    def __str__(self):
        s = str(self.N) + ','
        s += '(' + ','.join([str(p) for p in self.perms]) + '),'
        s += '(' + ','.join([str(o) for o in self.origin]) + '),'
        s += '(' + ','.join([str(r) for r in self.refinement]) + '),'
        s += '(' + ','.join([str(c) for c in self.coarsening]) + '),'
        s += '(' + ','.join([str(e) for e in self.extents]) + ')'
        return '[' + s + ']'
    
class ParamBuilder(object):
    
    def __init__(self, N):
        self.N = N
        self.perms = None
        self.origin = None
        self.refinement = None
        self.coarsening = None
        self.extents = None

    def to_params(self):
        if self.N != len(self.perms):
            print(f'Got length {len(self.perms)} for perms. Expected {self.N}.')
            exit(1)
        if self.N != len(self.origin):
            print(f'Got length {len(self.origin)} for origin. Expected {self.N}.')
            exit(1)
        if self.N != len(self.refinement):
            print(f'Got length {len(self.refinement)} for refinement. Expected {self.N}.')
            exit(1)
        if self.N != len(self.coarsening):
            print(f'Got length {len(self.coarsening)} for coarsening. Expected {self.N}.')
            exit(1)
        if self.N != len(self.extents):
            print(f'Got length {len(self.extents)} for extents. Expected {self.N}.')
            exit(1)
        return Params(self.N, self.perms, self.origin, self.refinement, self.coarsening, self.extents)

    def with_N(self, N):
        self.N = N
        return N

    def with_perms(self, perms):
        self.perms = perms
        return self

    def with_origin(self, origin):
        self.origin = origin
        return self

    def with_refinement(self, refinement):
        self.refinement = refinement        
        return self

    def with_coarsening(self, coarsening):
        self.coarsening = coarsening
        return self

    def with_extents(self, extents):
        self.extents = extents
        return self

    @staticmethod
    def defaults(N):
        perms = []
        for i in range(N):
            perms.append(i)
        origin = [0] * N
        refinement = [1] * N
        coarsening = [1] * N
        extents = [1] * N
        pbuild = ParamBuilder(N)
        pbuild.with_perms(perms)
        pbuild.with_origin(origin)
        pbuild.with_refinement(refinement)
        pbuild.with_coarsening(coarsening)
        pbuild.with_extents(extents)
        return pbuild


class TSpace(object):
    
    def __init__(self):
        self.is_universe = False
        self.is_block = False
        self.params = None
        self.parent = None
        self.children = set()

    def __init__(self, is_universe, is_block, params, parent):
        self.is_universe = is_universe
        self.is_block = is_block
        self.params = params
        self.parent = parent
        self.children = set()
        if not is_universe:
            if parent is None:
                print('Non-universe space must have a non-null parent!')
                exit(1)
            self.parent.children.add(self)

    def __str__(self):
        s = f'{self.is_universe},{self.is_block},{self.params}'
        return s

    def _get_upstream_path(self, parents):
        if self.parent is not None:
            parents.append(self.parent)
            self.parent._get_upstream_path(parents)

    def get_upstream_path(self):
        parents = []
        self._get_upstream_path(parents)
        return parents

    def Phi(self, c_L, F, L, R):
        if L >= R:
            c_r = []
            for i in range(L-R,L):
                c_r.append(c_L[i])
            return c_r
        else:
            c_r = [0] * (R-L)
            for i in range(L):
                c_r.append(c_L[i])
            return c_r

    # move one
    def phi_upstream_mapping(self, c_z):
        if len(c_z) != self.params.N:
            print(f"Expected c_z to have {self.params.N} components. Got {len(c_z)}")
            exit(1)
        if self.parent is None:
            if not self.is_universe:
                print("Only the universe can exist without a parent!")
                exit(1)
            return c_z
        else:
            rho = lambda t, k : int(((t+self.params.origin[k])*self.params.coarsening[k])/self.params.refinement[k])
            c_p = [rho(c_k,k) for k,c_k in enumerate(c_z)]
            # permute
            c_pp = []
            for i in range(self.params.N):
                for k in range(self.params.N):
                    if self.params.perms[k] == i:
                        c_pp.append(c_p[k])
                        break
            c_y = self.Phi(c_pp,0,self.params.N,self.parent.params.N)
            return c_y

    # move one
    def phi_downstream_mapping(self, c_y):
        if self.is_universe == True:
            print("Cannot perform a downstream mapping on the universe!")
            exit(1)
        if len(c_y) != self.parent.params.N:
            print(f"Expected c_y to have {self.parent.params.N} components. Got {len(c_y)}")
        rho = lambda t, k : int((t*self.params.refinement[k]/self.params.coarsening[k]))-self.params.origin[k]
        c_p = self.Phi(c_y,0,self.parent.params.N,self.params.N)
        c_pp = []
        for i in range(self.params.N):
            c_pp.append(c_p[self.params.perms[i]])
        c_z = [rho(c_k,k) for k,c_k in enumerate(c_pp)]
        return c_z        

    def access(self, c_z):
        if len(c_z) != self.params.N:
            print(f'Expected c_z to have {self.params.N} components. Got {len(c_z)}')
            exit(1)
        if self.is_block or self.is_universe:
            return c_z
        else:
            if self.parent is None:
                printf('Expected parent')
                exit(1)
            return self.parent.access(self.phi_upstream_mapping(c_z))

    def block_copy(self):
        params = ParamBuilder.defaults(self.params.N)
        params = params.with_extents(self.params.extents).to_params()
        child = TSpace(False, True, params, self)
        self.children.add(child)
        return child

    def view_copy(self):
        params = ParamBuilder.defaults(self.params.N)
        params = params.with_extents(self.params.extents).to_params()
        child = TSpace(False, False, params, self)
        self.children.add(child)
        return child
    
    def partition(self, *pparams):
        the_params = [p for p in pparams]
        if len(the_params) != self.params.N:
            print(f'Expected partition to have {self.params.N} parameters. Got {len(the_params)}')
            exit(1)
        for p in the_params:
            if len(p) != 3:
                print(f'Expected 3 values per partition parameter. Got {len(p)}')
                exit(1)
        V_p_z = TSpace(False, False, ParamBuilder.defaults(self.params.N).with_origin([p[0] for p in the_params]).to_params(), self)
        self.children.add(V_p_z)
        V_pp_z = TSpace(False, False, ParamBuilder.defaults(self.params.N).with_coarsening([p[2] for p in the_params]).to_params(), V_p_z)
        V_p_z.children.add(V_pp_z)
        V_z = TSpace(False, False, ParamBuilder.defaults(self.params.N).with_extents([int((p[1]-p[0]-1)/p[2])+1 for p in the_params]).to_params(), V_pp_z)
        V_pp_z.children.add(V_z)
        return V_z

    def slice(self, K):
        if K <= 0 or K > self.params.N:
            print(f'Expected K to be 0 < K <= {self.params.N}. Got {K}')
            exit(1)
        child = TSpace(False, False, ParamBuilder.defaults(K).with_extents(self.Phi(self.params.extents, 0, self.params.N, K)).to_params(), self)
        self.children.add(child)
        return child

    def vpermute(self, F):
        if len(F) != self.params.N:
            print(f'Expected {self.params.N} permutation factors. Got {len(F)}.')
            exit(1)
        found = set()
        for f in F:
            if f < 0:
                print(f'Invalid permutation factor {f}. Must be >= 0')
                exit(1)
            if f in found:
                print(f'Duplicate permutation factor {f}.')
                exit(1)
            found.add(f)
        E_v = []
        for f in F:
            E_v.append(self.params.extents[f])
        child = TSpace(False, False, ParamBuilder.defaults(self.params.N).with_extents(E_v).with_perms(F).to_params(), self)
        self.children.add(child)
        return child

    def vrefine(self, F):
        if len(F) != self.params.N:
            print(f'Expected {self.params.N} refinement factors. Got {len(F)}.')
            exit(1)
        for f in F:
            if f < 0:
                print(f'Invalid refinement factor {f}. Must be >= 0')
                exit(1)
        E_v = [self.params.extents[i] * F[i] for i in range(self.params.N)]
        child = TSpace(False, False, ParamBuilder.defaults(self.params.N).with_extents(E_v).with_refinement(F).to_params(), self)        
        self.children.add(child)
        return child

    def vcoarsen(self, F):
        if len(F) != self.params.N:
            print(f'Expected {self.params.N} refinement factors. Got {len(F)}.')
            exit(1)
        for f in F:
            if f < 0:
                print(f'Invalid refinement factor {f}. Must be >= 0')
                exit(1)
        E_v = [int((self.params.extents[i] - 1)/ F[i])+1 for i in range(self.params.N)]
        child = TSpace(False, False, ParamBuilder.defaults(self.params.N).with_extents(E_v).with_coarsening(F).to_params(), self)        
        self.children.add(child)
        return child

# Manual builds
universe = TSpace(True, False, ParamBuilder.defaults(3).to_params(), None)
block_A = TSpace(False, True, ParamBuilder.defaults(3).with_extents([2,4,10]).to_params(), universe)
block_B = TSpace(False, True, ParamBuilder.defaults(2).with_extents([6,8]).to_params(), universe)
block_C = TSpace(False, True, ParamBuilder.defaults(2).with_extents([6,8]).with_origin([1,2]).to_params(), universe)
block_D = TSpace(False, True, ParamBuilder.defaults(3).with_extents([2,4,10]).with_perms([2,0,1]).to_params(), universe)
block_E = TSpace(False, True, ParamBuilder.defaults(3).with_extents([1,2,10]).with_coarsening([2,2,1]).to_params(), universe)
block_F = TSpace(False, True, ParamBuilder.defaults(3).with_extents([4,16,10]).with_refinement([2,3,1]).to_params(), universe)
# Using creators
block_G = block_F.block_copy()
view_H = block_F.view_copy()
view_I = block_C.partition([0,4,2],[1,3,1])
view_J = block_E.slice(2)
view_K = block_E.slice(1)

print(universe)
print(block_A)
print(block_B)
print(block_C)
print(block_D)
print(block_E)
print(view_I)
print(view_J)
print(view_K)

#print("Moving point (0,0,0) from block_A to universe")
#p = block_A.phi_upstream_mapping([0,0,0])
#assert p == [0,0,0]
#
#print("Moving point (1,2,3) from block_A to universe")
#p = block_A.phi_upstream_mapping([1,2,3])
#assert p == [1,2,3]
#
#print("Moving point (0,0) from block_B to universe")
#p = block_B.phi_upstream_mapping([0,0])
#assert p == [0,0,0]
#
#print("Moving point (1,2) from block_B to universe")
#p = block_B.phi_upstream_mapping([1,2])
#assert p == [0,1,2]
#
#print("Moving point (0,0) from block_C to universe")
#p = block_C.phi_upstream_mapping([0,0])
#assert p == [0,1,2]
#
#print("Moving point (-2,4) from block_C to universe")
#p = block_C.phi_upstream_mapping([-2,4])
#assert p == [0,-1,6]
#
#print("Moving point (3,1,2) from block_D to universe")
#p = block_D.phi_upstream_mapping([3,1,2])
#assert p == [1,2,3]
#
#print("Moving point (1,2,3) from block_E to universe")
#p = block_E.phi_upstream_mapping([1,2,3])
#assert p == [2,4,3]
#
#print("Moving point (1,3,4) from block_F to universe")
#p = block_F.phi_upstream_mapping([1,3,4])
#assert p == [0,1,4]
#
#print("Moving point (0,0,0) from universe to block_A")
#p = block_A.phi_downstream_mapping([0,0,0])
#assert p == [0,0,0]
#
#print("Moving point (0,0,0) from universe to block_B")
#p = block_B.phi_downstream_mapping([0,0,0])
#assert p == [0,0]
#
#print("Moving point (0,0,0) from universe to block_C")
#p = block_C.phi_downstream_mapping([0,0,0])
#assert p == [-1,-2]
#
#print("Moving point (17,1,2) from universe to block_C")
#p = block_C.phi_downstream_mapping([17,1,2])
#assert p == [0,0]
#
#print("Moving point (1,2,3) from universe to block_D")
#p = block_D.phi_downstream_mapping([1,2,3])
#assert p == [3,1,2]
#
#print("Moving point (4,5,3) from universe to block_E")
#p = block_E.phi_downstream_mapping([4,5,3])
#assert p == [2,2,3]
#
#print("Moving point (4,5,3) from universe to block_F")
#p = block_F.phi_downstream_mapping([4,5,3])
#assert p == [8,15,3]
#
#print("Accessing point (3,4,5) in block_F")
#p = block_F.access([3,4,5])
#assert p == [3,4,5]
#
#print("Moving point (4,5,3) from block_F to block_G")
#p = block_G.phi_downstream_mapping([4,5,3])
#assert p == [4,5,3]
#
#print("Moving point (4,5,3) from block_F to view_H")
#p = view_H.phi_downstream_mapping([4,5,3])
#assert p == [4,5,3]
#
#print("Accessing point (0,0) in view_I")
#p = view_I.access([0,0])
#assert p == [0,1]
#
#print("Accessing point (0,1) in view_I")
#p = view_I.access([0,1])
#assert p == [0,2]
#
#print("Accessing point (1,0) in view_I")
#p = view_I.access([1,0])
#assert p == [2,1]
#
#print("Accessing point (1,1) in view_I")
#p = view_I.access([1,1])
#assert p == [2,2]
#
#print("Accessing point (3,5) in view_J")
#p = view_J.access([3,5])
#assert p == [0,3,5]
#
#print("Accessing point (5,) in view_K")
#p = view_K.access([5])
#assert p == [0,0,5]
#
view_L = view_H.vpermute([1,2,0])
print(view_L)

print('Accessing point (0,1,2) in view_L')
p = view_L.access([0,1,2])
assert p == [2,0,1]

view_M = view_L.vrefine([2,3,4])
print(view_M)

print('Accessing point (0,3,8) in view_M')
p = view_M.access([0,3,8])
assert p == [2,0,1]

view_N = view_M.vcoarsen([1,1,2])
print(view_N)

print('Accessing point (0,3,4) in view_N')
p = view_N.access([0,3,4])
assert p == [2,0,1]

frame = TSpace(False, True, ParamBuilder.defaults(2).with_extents([32,32]).to_params(), universe)
vframe = frame.vcoarsen([16,16])
assert vframe.params.extents == [2,2]
metadata = vframe.block_copy()
assert vframe.params.extents == [2,2]
vmetadata = metadata.vrefine([16,16])
print(frame)
print(vframe)
print(metadata)
print(vmetadata)
for i in range(16):
    for j in range(16):
        assert(vmetadata.access([i,j]) == [0,0])
        assert(vmetadata.access([i+16,j]) == [1,0])
        assert(vmetadata.access([i,j+16]) == [0,1])
        assert(vmetadata.access([i+16,j+16]) == [1,1])
