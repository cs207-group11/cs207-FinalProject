import reaction_coeffs

def test_const_result():
    assert reaction_coeffs.const(10) == 10

def test_arr_value_A():
    try:
        reaction_coeffs.arr(-1,10**3,10**2)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_arr_value_T():
    try:
        reaction_coeffs.arr(10**7,10**3,-1)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_arr_value_R():
    try:
        reaction_coeffs.arr(1,0.5,10**3,-1)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_mod_arr_value_A():
    try:
        reaction_coeffs.mod_arr(-1,0.5,10**3,10**2)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_mod_arr_value_T():
    try:
        reaction_coeffs.mod_arr(10**7,0.5,10**3,-1)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_mod_arr_value_R():
    try:
        reaction_coeffs.mod_arr(1,0.5,10**3,10**2,-1)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_mod_arr_value_b():
    try:
        reaction_coeffs.mod_arr(10**7,"1",10**3,1)
    except TypeError as err:
        assert(type(err) == TypeError)