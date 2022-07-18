import dill

def serialize_fit(fit,filename):
    if '__dll_func__' in fit.__model__.__dict__:
        fit.__model__.__dll_func__ = None # remove lambda function to make Pickle work
    if '__initialized_func__' in fit.__model__.__dict__:
        fit.__model__.__initialized_func__ = False # remove __initialized_func__ to make Pickle work
    if '__initialized_jac__' in fit.__model__.__dict__:
        fit.__model__.__initialized_jac__ = False # remove __initialized_jac__ to make Pickle work
    if '__lambdified_func__' in fit.__model__.__dict__:
        fit.__model__.__lambdified_func__ = None # remove __lambdified_func__ to make Pickle work    
    if '__numbified_func__' in fit.__model__.__dict__:
        fit.__model__.__numbified_func__ = None # remove __numbified_func__ to make Pickle work    
    if '__lambdified_jac__' in fit.__model__.__dict__:
        fit.__model__.__lambdified_jac__ = None # remove __lambdified_jac__ to make Pickle work    
    if '__numbified_jac__' in fit.__model__.__dict__:
        fit.__model__.__numbified_jac__ = None # remove __numbified_jac__ to make Pickle work   
    if '__function_fortran__' in fit.__model__.__dict__:
        fit.__model__.__function_fortran__ = None # remove __function_fortran__ to make Pickle work           
    with open(filename, 'wb') as output:
        dill.dump(fit, output, dill.DEFAULT_PROTOCOL)
       
def serialize_model(model,filename):
    with open(filename,'wb') as fil:   
        if '__dll_func__' in model.__dict__:
            model.__dll_func__ = None # remove lambda function to make Pickle work
        if '__initialized_func__' in model.__dict__:
            model.__initialized_func__ = False # remove __initialized_func__ to make Pickle work
        if '__initialized_jac__' in model.__dict__:
            model.__initialized_jac__ = False # remove __initialized_jac__ to make Pickle work
        if '__lambdified_func__' in model.__dict__:
            model.__lambdified_func__ = None # remove __lambdified_func__ to make Pickle work    
        if '__numbified_func__' in model.__dict__:
            model.__numbified_func__ = None # remove __numbified_func__ to make Pickle work    
        if '__lambdified_jac__' in model.__dict__:
            model.__lambdified_jac__ = None # remove __lambdified_jac__ to make Pickle work    
        if '__numbified_jac__' in model.__dict__:
            model.__numbified_jac__ = None # remove __numbified_jac__ to make Pickle work    
        dill.dump(model,fil,dill.DEFAULT_PROTOCOL)

def deserialize_fit(filename):
    with open(filename, 'rb') as input:
        return dill.load(input)

def deserialize_model(filename):
    with open(filename, 'rb') as input:
        return dill.load(input)
