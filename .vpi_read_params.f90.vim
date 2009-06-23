let b:file_mod_time="1171583725"
syn keyword FortranFunction vpi_read_params
syn keyword FortranFunction vpi_init_sim_control
let b:num_functions = 2
let b:string = ['vpi_read_params', 'vpi_init_sim_control:vpi_read_params']
let b:type = ['program', 'subroutine']
let b:start_line = [1, 13]
let b:end_line = [23, 21]
let b:start_line_bracket = [2, 14]
let b:direct_modules_used_list = []
let b:direct_modules_use_all = []
let b:direct_modules_only_dict = {}
let b:direct_rename_dict = {}
let b:complete_module_dict = {}
syn match FortranFunctionArgument display /\%>12l\%<22l%\@<!\<pname\>/
syn match FortranFunctionArgument display /\%>12l\%<22l%\@<!\<pval\>/