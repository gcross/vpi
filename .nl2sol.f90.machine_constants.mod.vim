let b:module_direct_dict['machine_constants'] = {'interface:imdcon': {'arg:order': ['k'], 'k': ['INTEGER', 'INTENT(in)']}, 'imdcon': ['INTEGER', 'function', 'module:machine_constants', 'nl2sol.f90'], 'dp': ['INTEGER', 'PARAMETER', '=selected_real_kind(precision(1.d0),range(1.d0))', 'variable', 'module:machine_constants', 'nl2sol.f90'], 'interface:errstop': {'routine': ['CHARACTER(*)', 'INTENT(in)'], 'arg:order': ['routine', 'error'], 'error': ['CHARACTER(*)', 'INTENT(in)']}, 'errstop': ['subroutine', 'module:machine_constants', 'nl2sol.f90'], 'machine_constants': ['module', 'nl2sol.f90'], 'interface:rmdcon': {'arg:order': ['k'], 'k': ['INTEGER', 'INTENT(in)']}, 'rmdcon': ['REAL(dp)', 'function', 'module:machine_constants', 'nl2sol.f90']}
let b:default_accessibility['machine_constants'] = 'public'
let b:modules_used_list['machine_constants'] = []
let b:modules_use_all['machine_constants'] = []
let b:modules_only_dict['machine_constants'] = {}
let b:rename_dict['machine_constants'] = {}
let b:indirect_module_dict['machine_constants'] = {}