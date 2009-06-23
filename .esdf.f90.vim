let b:file_mod_time="1202406289"
syn keyword FortranFunction paramdf_init
syn keyword FortranFunction paramdf_string
syn keyword FortranFunction paramdf_integer
syn keyword FortranFunction paramdf_single
syn keyword FortranFunction paramdf_double
syn keyword FortranFunction paramdf_physical
syn keyword FortranFunction paramdf_defined
syn keyword FortranFunction paramdf_boolean
syn keyword FortranFunction paramdf_block
syn keyword FortranFunction paramdf_reduce
syn keyword FortranFunction paramdf_convfac
syn keyword FortranFunction paramdf_unit
syn keyword FortranFunction paramdf_file
syn keyword FortranFunction paramdf_newfile
syn keyword FortranFunction paramdf_lblchk
syn keyword FortranFunction paramdf_help
syn keyword FortranFunction paramdf_die
syn keyword FortranFunction paramdf_warn
syn keyword FortranFunction paramdf_warnout
syn keyword FortranFunction paramdf_close
syn keyword FortranFunction paramdf_dump
syn keyword FortranFunction help_system
let b:num_functions = 22
let b:string = ['paramdf_init', 'paramdf_string', 'paramdf_integer', 'paramdf_single', 'paramdf_double', 'paramdf_physical', 'paramdf_defined', 'paramdf_boolean', 'paramdf_block', 'paramdf_reduce', 'paramdf_convfac', 'paramdf_unit', 'paramdf_file', 'paramdf_newfile', 'paramdf_lblchk', 'paramdf_help', 'paramdf_die', 'paramdf_warn', 'paramdf_warnout', 'paramdf_close', 'paramdf_dump', 'help_system']
let b:type = ['subroutine', 'function', 'function', 'function', 'function', 'function', 'function', 'function', 'function', 'function', 'function', 'function', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine']
let b:start_line = [230, 375, 412, 452, 492, 533, 578, 618, 673, 727, 772, 811, 833, 853, 870, 899, 1030, 1039, 1049, 1074, 1082, 1153]
let b:end_line = [372, 409, 449, 489, 530, 575, 615, 670, 724, 769, 808, 830, 850, 867, 896, 1027, 1036, 1046, 1071, 1079, 1150, 1184]
let b:start_line_bracket = [231, 376, 413, 453, 493, 534, 579, 619, 674, 728, 773, 812, 834, 854, 871, 900, 1031, 1040, 1050, 1075, 1083, 1154]
syn keyword FortranModule paramdf
syn match FortranModuleVariable display /%\@<!\<llength\>/
syn match FortranModuleVariable display /%\@<!\<nphys\>/
syn match FortranModuleVariable display /%\@<!\<ndump\>/
syn match FortranModuleVariable display /%\@<!\<nrecords\>/
syn match FortranModuleVariable display /%\@<!\<nwarns\>/
syn match FortranModuleVariable display /%\@<!\<ndmp\>/
syn match FortranModuleVariable display /%\@<!\<keywords_not_in_list\>/
syn match FortranModuleVariable display /%\@<!\<llist\>/
syn match FortranModuleVariable display /%\@<!\<warns\>/
syn match FortranModuleVariable display /%\@<!\<dump\>/
syn match FortranModuleVariable display /%\@<!\<tlist\>/
syn match FortranModuleVariable display /%\@<!\<block_data\>/
syn match FortranModuleVariable display /%\@<!\<phy\>\%(([^)]*)\)*%\a\w*\>\%(([^)]*)\)*%\a\w*\>\|%\@<!\<phy\>\%(([^)]*)\)*%\a\w*\>\|%\@<!\<phy\>/
let b:direct_modules_used_list = ['dsp', 'paramdf_key', 'utilities', 'store']
let b:direct_modules_use_all = ['dsp', 'paramdf_key', 'utilities']
let b:direct_modules_only_dict = {'store': ['o']}
let b:direct_rename_dict = {}
let b:complete_module_dict = {}
syn keyword FortranTypeDef phys_unit
syn match FortranFunctionArgument display /\%>229l\%<373l%\@<!\<filename\>/
syn match FortranFunctionArgument display /\%>374l\%<410l%\@<!\<label\>/
syn match FortranFunctionArgument display /\%>374l\%<410l%\@<!\<default\>/
syn match FortranFunctionArgument display /\%>411l\%<450l%\@<!\<label\>/
syn match FortranFunctionArgument display /\%>411l\%<450l%\@<!\<default\>/
syn match FortranFunctionArgument display /\%>451l\%<490l%\@<!\<label\>/
syn match FortranFunctionArgument display /\%>451l\%<490l%\@<!\<default\>/
syn match FortranFunctionArgument display /\%>491l\%<531l%\@<!\<label\>/
syn match FortranFunctionArgument display /\%>491l\%<531l%\@<!\<default\>/
syn match FortranFunctionArgument display /\%>532l\%<576l%\@<!\<label\>/
syn match FortranFunctionArgument display /\%>532l\%<576l%\@<!\<default\>/
syn match FortranFunctionArgument display /\%>532l\%<576l%\@<!\<dunit\>/
syn match FortranFunctionArgument display /\%>577l\%<616l%\@<!\<label\>/
syn match FortranFunctionArgument display /\%>577l\%<616l%\@<!\<type\>/
syn match FortranFunctionArgument display /\%>617l\%<671l%\@<!\<label\>/
syn match FortranFunctionArgument display /\%>617l\%<671l%\@<!\<default\>/
syn match FortranFunctionArgument display /\%>672l\%<725l%\@<!\<label\>/
syn match FortranFunctionArgument display /\%>672l\%<725l%\@<!\<nlines\>/
syn match FortranFunctionArgument display /\%>726l\%<770l%\@<!\<string_untrimmed\>/
syn match FortranFunctionArgument display /\%>771l\%<809l%\@<!\<from\>/
syn match FortranFunctionArgument display /\%>771l\%<809l%\@<!\<to\>/
syn match FortranFunctionArgument display /\%>810l\%<831l%\@<!\<ierr\>/
syn match FortranFunctionArgument display /\%>832l\%<851l%\@<!\<unit\>/
syn match FortranFunctionArgument display /\%>832l\%<851l%\@<!\<filename\>/
syn match FortranFunctionArgument display /\%>832l\%<851l%\@<!\<ierr\>/
syn match FortranFunctionArgument display /\%>852l\%<868l%\@<!\<unit\>/
syn match FortranFunctionArgument display /\%>852l\%<868l%\@<!\<filename\>/
syn match FortranFunctionArgument display /\%>852l\%<868l%\@<!\<ierr\>/
syn match FortranFunctionArgument display /\%>869l\%<897l%\@<!\<string\>/
syn match FortranFunctionArgument display /\%>869l\%<897l%\@<!\<typ\>/
syn match FortranFunctionArgument display /\%>898l\%<1028l%\@<!\<helpword_in\>/
syn match FortranFunctionArgument display /\%>898l\%<1028l%\@<!\<searchword_in\>/
syn match FortranFunctionArgument display /\%>1029l\%<1037l%\@<!\<string\>/
syn match FortranFunctionArgument display /\%>1038l\%<1047l%\@<!\<string\>/
syn match FortranFunctionArgument display /\%>1081l\%<1151l%\@<!\<filename\>/
syn match FortranFunctionArgument display /\%>1152l\%<1185l%\@<!\<hword\>/
syn match FortranFunctionArgument display /\%>1152l\%<1185l%\@<!\<sword\>/