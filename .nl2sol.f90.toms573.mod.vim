let b:module_direct_dict['toms573'] = {'interface:livmul': {'arg:order': ['n', 'x', 'l', 'y'], 'x': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'y': ['REAL(dp)', 'INTENT(in)', '(:)'], 'l': ['REAL(dp)', 'INTENT(in)', '(:)'], 'n': ['INTEGER', 'INTENT(in)']}, 'livmul': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'nl2sol': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:ltsqar': {'a': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'arg:order': ['n', 'a', 'l'], 'l': ['REAL(dp)', 'INTENT(in)', '(:)'], 'n': ['INTEGER', 'INTENT(in)']}, 'ltsqar': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:vcopy': {'p': ['INTEGER', 'INTENT(in)'], 'arg:order': ['p', 'y', 'x'], 'x': ['REAL(dp)', 'INTENT(in)', '(:)'], 'y': ['REAL(dp)', 'INTENT(inout)', '(:)']}, 'interface:dotprd': {'p': ['INTEGER', 'INTENT(in)'], 'arg:order': ['p', 'x', 'y'], 'x': ['REAL(dp)', 'INTENT(in)', '(:)'], 'y': ['REAL(dp)', 'INTENT(in)', '(:)']}, 'dotprd': ['REAL(dp)', 'function', 'module:toms573', 'nl2sol.f90'], 'toms573': ['module', 'nl2sol.f90'], 'interface:assess': {'p': ['INTEGER', 'INTENT(in)'], 'x0': ['REAL(dp)', 'INTENT(inout)'], 'iv': ['INTEGER', 'INTENT(inout)', '(:)'], 'd': ['REAL(dp)', 'INTENT(in)', '(:)'], 'stlstg': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'v': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'arg:order': ['d', 'iv', 'p', 'step', 'stlstg', 'v', 'x', 'x0'], 'x': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'step': ['REAL(dp)', 'INTENT(inout)', '(:)']}, 'assess': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:itsmry': {'p': ['INTEGER', 'INTENT(in)'], 'iv': ['INTEGER', 'INTENT(inout)', '(:)'], 'd': ['REAL(dp)', 'INTENT(in)', '(:)'], 'v': ['INTEGER', 'INTENT(in)'], 'x': ['INTEGER', 'INTENT(in)'], 'arg:order': ['d', 'iv', 'p', 'v', 'x']}, 'itsmry': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:v2norm': {'p': ['INTEGER', 'INTENT(in)'], 'arg:order': ['p', 'x'], 'x': ['REAL(dp)', 'INTENT(in)', '(:)']}, 'slupdt': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:reldst': {'p': ['INTEGER', 'INTENT(in)'], 'x0': ['REAL(dp)', 'INTENT(in)', '(:)'], 'd': ['REAL(dp)', 'INTENT(in)', '(:)'], 'x': ['REAL(dp)', 'INTENT(in)', '(:)'], 'arg:order': ['p', 'd', 'x', 'x0']}, 'reldst': ['REAL(dp)', 'function', 'module:toms573', 'nl2sol.f90'], 'interface:qapply': {'p': ['INTEGER', 'INTENT(in)'], 'r': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'arg:order': ['n', 'p', 'j', 'r', 'ierr'], 'ierr': ['REAL(dp)', 'INTENT(inout)'], 'j': ['REAL(dp)', 'INTENT(in)', '(:,:)'], 'n': ['INTEGER', 'INTENT(in)']}, 'qapply': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:nl2sol': {'p': ['INTEGER', 'INTENT(in)'], 'uiparm': ['INTEGER', 'INTENT(inout)', 'OPTIONAL', '(:)'], 'iv': ['INTEGER', 'INTENT(inout)', '(:)'], 'arg:order': ['n', 'p', 'x', 'iv', 'v', 'uiparm', 'urparm', 'ufparm'], 'v': ['INTEGER', 'INTENT(inout)'], 'x': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'urparm': ['REAL(dp)', 'INTENT(inout)', 'OPTIONAL', '(:)'], 'n': ['INTEGER', 'INTENT(in)'], 'ufparm': ['REAL(dp)', 'INTENT(inout)', 'OPTIONAL']}, 'interface:litvmu': {'arg:order': ['n', 'x', 'l', 'y'], 'x': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'y': ['REAL(dp)', 'INTENT(in)', '(:)'], 'l': ['REAL(dp)', 'INTENT(in)', '(:)'], 'n': ['INTEGER', 'INTENT(in)']}, 'litvmu': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:stopx': {'arg:order': []}, 'interface:covclc': {'p': ['INTEGER', 'INTENT(in)'], 'r': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'iv': ['REAL(dp)', 'INTENT(in)'], 'd': ['REAL(dp)', 'INTENT(in)', '(:)'], 'v': ['REAL(dp)', 'INTENT(inout)'], 'x': ['REAL(dp)', 'INTENT(inout)'], 'arg:order': ['covirc', 'd', 'iv', 'j', 'n', 'p', 'r', 'v', 'x'], 'covirc': ['INTEGER', 'INTENT(inout)'], 'j': ['REAL(dp)', 'INTENT(inout)', '(:,:)'], 'n': ['INTEGER', 'INTENT(in)']}, 'interface:vscopy': {'p': ['INTEGER', 'INTENT(in)'], 's': ['REAL(dp)', 'INTENT(in)'], 'arg:order': ['p', 'y', 's'], 'y': ['REAL(dp)', 'INTENT(inout)', '(:)']}, 'vscopy': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'lsqrt': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:slvmul': {'p': ['INTEGER', 'INTENT(in)'], 's': ['REAL(dp)', 'INTENT(in)', '(:)'], 'arg:order': ['p', 'y', 's', 'x'], 'x': ['REAL(dp)', 'INTENT(in)', '(:)'], 'y': ['REAL(dp)', 'INTENT(inout)', '(:)']}, 'slvmul': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:dupdat': {'p': ['REAL(dp)', 'INTENT(in)'], 'iv': ['INTEGER', 'INTENT(in)', '(:)'], 'd': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'v': ['REAL(dp)', 'INTENT(in)', '(:)'], 'arg:order': ['d', 'iv', 'j', 'n', 'p', 'v'], 'j': ['REAL(dp)', 'INTENT(in)', '(:,:)'], 'n': ['REAL(dp)', 'INTENT(in)']}, 'dupdat': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:vaxpy': {'p': ['INTEGER', 'INTENT(in)'], 'a': ['REAL(dp)', 'INTENT(in)'], 'arg:order': ['p', 'w', 'a', 'x', 'y'], 'w': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'x': ['REAL(dp)', 'INTENT(in)', '(:)'], 'y': ['REAL(dp)', 'INTENT(in)', '(:)']}, 'vcopy': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:nl2itr': {'p': ['INTEGER', 'INTENT(in)'], 'jac': ['INTEGER', 'INTENT(inout)'], 'iv': ['INTEGER', 'INTENT(inout)', '(:)'], 'd': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'v': ['INTEGER', 'INTENT(in)'], 'arg:order': ['d', 'iv', 'jac', 'n', 'nn', 'p', 'r', 'v', 'x'], 'x': ['INTEGER', 'INTENT(in)'], 'nn': ['INTEGER', 'INTENT(in)'], 'r': ['INTEGER', 'INTENT(in)'], 'n': ['INTEGER', 'INTENT(in)']}, 'nl2itr': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:lmstep': {'p': ['INTEGER', 'INTENT(in)'], 'r': ['INTEGER', 'INTENT(in)'], 'd': ['REAL(dp)', 'INTENT(in)', '(:)'], 'w': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'v': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'g': ['REAL(dp)', 'INTENT(in)', '(:)'], 'ka': ['INTEGER', 'INTENT(inout)'], 'step': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'ipivot': ['INTEGER', 'INTENT(inout)', '(:)'], 'arg:order': ['d', 'g', 'ierr', 'ipivot', 'ka', 'p', 'qtr', 'r', 'step', 'v', 'w'], 'ierr': ['INTEGER', 'INTENT(inout)'], 'qtr': ['INTEGER', 'INTENT(in)']}, 'lmstep': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:linvrt': {'arg:order': ['n', 'lin', 'l'], 'lin': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'l': ['REAL(dp)', 'INTENT(in)', '(:)'], 'n': ['INTEGER', 'INTENT(in)']}, 'linvrt': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:gqtstp': {'p': ['INTEGER', 'INTENT(in)'], 'ka': ['INTEGER', 'INTENT(inout)'], 'd': ['REAL(dp)', 'INTENT(in)', '(:)'], 'v': ['INTEGER', 'INTENT(in)'], 'w': ['INTEGER', 'INTENT(in)'], 'dig': ['REAL(dp)', 'INTENT(in)', '(:)'], 'step': ['INTEGER', 'INTENT(in)'], 'l': ['INTEGER', 'INTENT(inout)'], 'dihdi': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'arg:order': ['d', 'dig', 'dihdi', 'ka', 'l', 'p', 'step', 'v', 'w']}, 'gqtstp': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:nl2sno': {'p': ['INTEGER', 'INTENT(in)'], 'uiparm': ['INTEGER', 'INTENT(inout)', 'OPTIONAL', '(:)'], 'iv': ['INTEGER', 'INTENT(inout)', '(:)'], 'arg:order': ['n', 'p', 'x', 'iv', 'v', 'uiparm', 'urparm', 'ufparm'], 'v': ['INTEGER', 'INTENT(inout)'], 'x': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'urparm': ['REAL(dp)', 'INTENT(inout)', 'OPTIONAL', '(:)'], 'n': ['INTEGER', 'INTENT(in)'], 'ufparm': ['REAL(dp)', 'INTENT(inout)', 'OPTIONAL']}, 'nl2sno': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:slupdt': {'p': ['INTEGER', 'INTENT(in)'], 'a': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'step': ['INTEGER', 'INTENT(in)'], 'cosmin': ['REAL(dp)', 'INTENT(in)'], 'arg:order': ['a', 'cosmin', 'p', 'size', 'step', 'u', 'w', 'wchmtd', 'wscale', 'y'], 'w': ['INTEGER', 'INTENT(in)'], 'y': ['REAL(dp)', 'INTENT(out)'], 'wchmtd': ['INTEGER', 'INTENT(in)'], 'u': ['INTEGER', 'INTENT(in)'], 'wscale': ['REAL(dp)', 'INTENT(out)'], 'size': ['INTEGER', 'INTENT(in)']}, 'v2norm': ['REAL(dp)', 'function', 'module:toms573', 'nl2sol.f90'], 'vaxpy': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:lsvmin': {'p': ['INTEGER', 'INTENT(in)'], 'arg:order': ['p', 'l', 'x', 'y', 'fn_val'], 'fn_val': ['REAL(dp)', 'INTENT(out)'], 'x': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'y': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'l': ['REAL(dp)', 'INTENT(in)', '(:)']}, 'lsvmin': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'stop_nl2sol': ['LOGICAL', '=.false.', 'variable', 'module:toms573', 'nl2sol.f90'], 'interface:parchk': {'p': ['INTEGER', 'INTENT(in)'], 'iv': ['INTEGER', 'INTENT(inout)', '(:)'], 'nn': ['INTEGER', 'INTENT(in)'], 'v': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'arg:order': ['iv', 'n', 'nn', 'p', 'v'], 'n': ['INTEGER', 'INTENT(in)']}, 'parchk': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:qrfact': {'sum': ['INTEGER', 'INTENT(inout)'], 'arg:order': ['m', 'n', 'qr', 'alpha', 'ipivot', 'ierr', 'nopivk', 'sum'], 'ierr': ['INTEGER', 'INTENT(inout)'], 'alpha': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'qr': ['REAL(dp)', 'INTENT(inout)', '(:,:)'], 'ipivot': ['INTEGER', 'INTENT(inout)', '(:)'], 'm': ['INTEGER', 'INTENT(in)'], 'n': ['INTEGER', 'INTENT(in)'], 'nopivk': ['INTEGER', 'INTENT(inout)']}, 'qrfact': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'stopx': ['LOGICAL', 'function', 'module:toms573', 'nl2sol.f90'], 'covclc': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:rptmul': {'p': ['REAL(dp)', 'INTENT(in)'], 'func': ['INTEGER', 'INTENT(in)'], 'arg:order': ['func', 'ipivot', 'j', 'p', 'rd', 'x', 'y', 'z'], 'ipivot': ['INTEGER', 'INTENT(in)', '(:)'], 'x': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'y': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'j': ['REAL(dp)', 'INTENT(in)', '(:,:)'], 'z': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'rd': ['REAL(dp)', 'INTENT(in)', '(:)']}, 'rptmul': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:dfault': {'iv': ['INTEGER', 'INTENT(inout)', '(:)'], 'arg:order': ['iv', 'v'], 'v': ['REAL(dp)', 'INTENT(inout)', '(:)']}, 'dfault': ['subroutine', 'module:toms573', 'nl2sol.f90'], 'interface:lsqrt': {'a': ['REAL(dp)', 'INTENT(in)', '(:)'], 'arg:order': ['n1', 'n', 'l', 'a', 'irc'], 'irc': ['INTEGER', 'INTENT(out)'], 'n1': ['INTEGER', 'INTENT(in)'], 'l': ['REAL(dp)', 'INTENT(inout)', '(:)'], 'n': ['INTEGER', 'INTENT(in)']}}
let b:default_accessibility['toms573'] = 'public'
let b:modules_used_list['toms573'] = ['machine_constants']
let b:modules_use_all['toms573'] = ['machine_constants']
let b:modules_only_dict['toms573'] = {}
let b:rename_dict['toms573'] = {}
let b:indirect_module_dict['toms573'] = {'interface:imdcon': {'arg:order': ['k'], 'k': ['INTEGER', 'INTENT(in)']}, 'imdcon': ['INTEGER', 'function', 'module:machine_constants', 'nl2sol.f90'], 'dp': ['INTEGER', 'PARAMETER', '=selected_real_kind(precision(1.d0),range(1.d0))', 'variable', 'module:machine_constants', 'nl2sol.f90'], 'interface:errstop': {'routine': ['CHARACTER(*)', 'INTENT(in)'], 'arg:order': ['routine', 'error'], 'error': ['CHARACTER(*)', 'INTENT(in)']}, 'errstop': ['subroutine', 'module:machine_constants', 'nl2sol.f90'], 'machine_constants': ['module', 'nl2sol.f90'], 'interface:rmdcon': {'arg:order': ['k'], 'k': ['INTEGER', 'INTENT(in)']}, 'rmdcon': ['REAL(dp)', 'function', 'module:machine_constants', 'nl2sol.f90']}