import itertools

# Load the file to be tested
load("affinea2.sage")
heart = [P1, P2, P3]

# We currently don't have code to access braid group elements by name. So let's temporarily build a dictionary to keep track of the generators.
generators = {}

generators['s_1'] = s1
generators['s_2'] = s2
generators['s_3'] = s3
generators['t_1'] = t1
generators['t_2'] = t2
generators['t_3'] = t3

def get_gen_by_sgn(i, s):
    """
    Return the name of sigma_i^s where s is a sign in {+1, -1}.
    Returns either 's_i' or 't_i'.
    """
    if s == 1:
        return 's_' + str(i)
    else:
        return 't_' + str(i)
    
# Currently we focus on testing homological linearity for the root-theoretically linear braids output by a separate piece of code.
def name_to_braid(bname):
    """
    Given a list of names of generators, return the function that is the action of the braid on the category.
    """
    return composeAll([generators[x] for x in bname])

# def apply_gen(g, obj):
#     """
#     Apply named generator g to object obj and return the result.
#     """

# def check_linearity_on(bname, obj):
#     """
#     Check whether a named braid is linear on an object obj. That is, whether it keeps it within the +/-1 t-slice.
#     """
#     g = bname.pop()

def vector_and_sign_to_braid(v,e):
    """
    Given a vector v of elements in {1,2,3} and a sign vector e of the same length, returns a list of 
    either 's_i' (for positive e_i) or 't_i' (for negative e_i).
    """
    return [get_gen_by_sgn(i,s) for (i,s) in zip(v,e)]

def is_homologically_linear(bname):
    """
    Given the name bname of a braid, check whether it is
    homologically linear. That is, whether it sends the heart to the
    [-1,1] t-slice. Do this one at a time for each object in heart.
    Return true if bname is homologically linear.

    Recall that heart = [P1, P2, P3]
    """
    b = name_to_braid(bname)
    for x in heart:
        bx = b(x)
        in_pm_1 = bx.minLevel() >= -1 and bx.maxLevel() <= 1
        if not in_pm_1:
            return False
    return True

def check_homological_linearity(v, e):
    """
    Given a vector v of elements in {1,2,3} and a sign vector e of
    the same length, let b be the braid that is the product of
    sigma_i^e_i.

    Return true if b is homologically linear, and false otherwise.

    """
    bname = vector_and_sign_to_braid(v,e)
    return is_homologically_linear(bname)


def is_cancelling_pair(s,t):
    """
    Given two letter names of the form 's_i' or 't_j', return true iff they form a cancelling pair.
    """
    return s[-1] == t[-1] and s[0] != t[0]

# def cancelling_partner(l):
#     """
#     Given a letter such as 's_i' or 't_j', return its cancelling partner.
#     """
#     x,n = l[0], l[-1]
#     if x == 's':
#         return 't_' + n
#     else:
#         returs 's_' + n


def is_cancellable_braid_move(bname):
    """
    Given a braid name, return True if a braid move at the beginning of the input can change a -++ word to a ++- word such that the - will cancel with the next +.
    """
    if len(bname) <= 3:
        return False
    else:
        a,b,c,d = bname[0],bname[1],bname[2],bname[3]
        if not is_cancelling_pair(a,c):
            return False
        if b == d:
            return True
        return False

def is_cancellable_si_ti(bname,place):
    """
    Given a braid name, return True if you can cancel a pair of si/ti from the place'th place in the braid name.
    Otherwise, return False.
    """
    if len(bname) <= place + 1:
        return False
    else:
        return is_cancelling_pair(bname[place], bname[place + 1])

def prune_si_ti(bname):
    """
    Successively cancel pairs of si/ti from the braid to reduce it. This code can surely be improved, but it works.
    """
    output = bname
    count = 0
    while count in range(len(output)):
        for i in range(len(output)):
            while(is_cancellable_si_ti(output, i)):
                output = output[:i] + output[i+2:]
        count = count + 1
    return output

def is_pm_or_mp(bname):
    """
    Return True if the braid name is of the form b+ times b- or b- times b+. Otherwise return False.
    """
    groups = [x[0] for x in itertools.groupby(bname, lambda y : y[0])]
    return len(groups) <= 2

load("linear_outputs_for_sage.sage")

pruned_braids = {}
for k in fixed_signs.keys():
    candidate_vectors,fixed_sign = sign_vectors[k], fixed_signs[k]
    bnames = [vector_and_sign_to_braid(v,fixed_sign) for v in candidate_vectors]
    bnames_pruned = [prune_si_ti(b) for b in bnames if not is_cancellable_si_ti(b, 0)]
    bnames_pruned = [b for b in bnames_pruned if not is_cancellable_braid_move(b)]
    bnames_pruned = [b for b in bnames_pruned if not is_pm_or_mp(b)]

    pruned_braids[k] = bnames_pruned
    
non_homologically_linear_braids = {}
for k in pruned_braids.keys():
    print("-------------------------")    
    print("Starting check for ",k)
    output = []
    for b in pruned_braids[k]:
        if not is_homologically_linear(b):
            output = [b] + output
    non_homologically_linear_braids[k] = output
    print("Finished check for %s", k)
    print("-------------------------")

