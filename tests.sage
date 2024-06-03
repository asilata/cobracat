load("zigzagalgebra_element.pyx")
load("zigzagalgebra-cartan.sage")
load("projective-zigzagmodules-cartan.sage")
load("complexes-cartan.sage")
load("braidactions-cartan.sage")

def run_zz_algebra_profile(N=1000):
    import random
    Z = ZigZagAlgebra("A2",QQ)
    B = Z.basis()
    [ random.choice(B)*random.choice(B)  for i in range(N)]
    

def run_zz_algera_internals_test():
    cartan_types = ["A2", "D5", CartanType(['A', 2, 1]), CartanType(['A',1,1])]
    for ct in cartan_types:
        basis = _zz_basis(CartanType(ct))
        names = [_zz_get_name(b) for b in basis]
        for b in basis:
            print(b, _zz_deg(b))
        table = [_zz_right_multiplication_table(basis, b) for b in basis]
        print(basis, names, table)
        
    return

def run_zig_zag_algebra_test():
    cartan_types = ["A2", "D5", CartanType(['A', 2, 1]), CartanType(['A',1,1])]
    for ct in cartan_types:
        print("Trying to construct zig-zag algebra of Cartan type {}.".format(ct))
        try:
            Z = ZigZagAlgebra(ct, QQ)
            print("Constructed {}.".format(Z))
            TestSuite(Z).run(verbose=True)
        except:
            print("FAILED to construct zigzag algebra of Cartan type {}".format(ct))

def run_zig_zag_module_test():
    cartan_types = ["A2", "D5", CartanType(['A', 2, 1])]
    for ct in cartan_types:
        print("Constructing zig-zag algebra of Cartan type {}".format(ct))
        Z = ZigZagAlgebra(ct, QQ)
        print(Z)
        proj_modules = [ProjectiveZigZagModule(Z, v) for v in Z.vertices]
        P1, P2 = proj_modules[0], proj_modules[1]
        print(P1.hom(P2))
        assert P1.is_invertible(P1.idempotent)
        assert P1.is_annihilated_by(P2.idempotent)
        assert P1.invert(P1.idempotent) == P1.idempotent
    return

def minimize_using_matrix_test(k=QQ):
    ct = "A3"
    Z = ZigZagAlgebra(ct, k)
    Z.inject_variables()
    P1 = ProjectiveZigZagModule(Z, 1)
    P2 = ProjectiveZigZagModule(Z, 2)
    P3 = ProjectiveZigZagModule(Z, 3)
    C = ProjectiveComplex(Z)
    C.add_object_at(0, P1)
    C.add_object_at(0, P2)
    C.add_object_at(1, P1)
    C.add_map_at(0,0,0,e1)
    C.add_map_at(0,1,0,a21)
    D, P = C.minimize_using_matrix(with_qis=True)
    assert D.min_index == 0
    assert D.max_index == 0
    assert D.objects[0][0] == P2
    assert is_chain_map(D,C,P)
    return C, D, P

def minimization_profile_setup(k=QQ, minimize=False):
    Z = ZigZagAlgebra("A2", k)
    p = {}
    for i in [1,2]:
        fpi = ProjectiveZigZagModule(Z, i)
        pi = ProjectiveComplex(Z)
        pi.add_object_at(0, fpi)
        p[i] = pi
    s = {i : (lambda C, i = i: sigma(Z, i, C, minimize=minimize)) for i in [1,2]}
    return p,s

def sigma_new_setup(k=QQ, minimize=False):
    Z = ZigZagAlgebra("A2", k)
    p = {}
    for i in [1,2]:
        fpi = ProjectiveZigZagModule(Z, i)
        pi = ProjectiveComplex(Z)
        pi.add_object_at(0, fpi)
        p[i] = pi
    s = {i : (lambda C, i = i: sigma_new(C, i, Z, minimize=minimize)) for i in [1,2]}
    return p,s

def run_projective_complex_test():
    return

def hom_complex_test():
    Z = ZigZagAlgebra("A2", k=QQ)
    p = {}
    for i in [1,2]:
        fpi = ProjectiveZigZagModule(Z, i)
        pi = ProjectiveComplex(Z)
        pi.add_object_at(0, fpi)
        p[i] = pi
    h1 = p[1].hom(p[2])
    assert h1.min_index == 0
    assert h1.max_index == 0
    assert len(h1.objects[0]) == 1
    assert h1.objects[0][0].grade() == -1

    h2 = p[1].hom(p[1])    
    assert h2.min_index == 0
    assert h2.max_index == 0
    assert len(h2.objects[0]) == 2
    assert set([x.grade() for x in h2.objects[0]]) == set([0,-2])
    
    return

def hom_complex_test_with_homs():
    p,s = minimization_profile_setup(minimize=True)
    X = s[1](s[1](p[2]))
    Y = p[2]
    hxy, mxy = X.hom(Y, with_homs=True)
    hyx, myx = Y.hom(X, with_homs=True)
    assert hxy.min_index == 0
    assert hxy.max_index == 2
    assert hyx.min_index == -2
    assert hyx.max_index == 0
    assert set(mxy.keys()) == set([2,0])
    assert set(myx.keys()) == set([-2, 0])
    return

def run_braid_actions_test():
    return

def run_all_tests():
    print("Running test for zig-zag algebras.")
    run_zig_zag_algebra_test()
    print("Running test for zig-zag modules.")    
    run_zig_zag_module_test()
    print("Running test for projective complexes.")        
    run_projective_complex_test()
    print("Running test for braid actions.")            
    run_braid_actions_test()
    return
