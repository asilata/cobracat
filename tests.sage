load("zigzagalgebra-cartan.sage")
load("projective-zigzagmodules-cartan.sage")

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
    cartan_types = ["A2", "D5", CartanType(['A', 2, 1]), CartanType(['A',1,1])]
    for ct in cartan_types:
        print("Constructing zig-zag algebra of Cartan type {}".format(ct))
        Z = ZigZagAlgebra(ct, QQ)
        print(Z)
        proj_modules = [ProjectiveZigZagModule(Z, v) for v in Z.vertices]
        P1, P2 = proj_modules[0], proj_modules[1]
        print(P1.hom(P2))
    return

def run_projective_complex_test():
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
