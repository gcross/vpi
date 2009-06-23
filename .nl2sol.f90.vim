let b:file_mod_time="1231266904"
syn keyword FortranFunction imdcon
syn keyword FortranFunction rmdcon
syn keyword FortranFunction errstop
syn keyword FortranFunction nl2sol
syn keyword FortranFunction nl2sno
syn keyword FortranFunction nl2itr
syn keyword FortranFunction assess
syn keyword FortranFunction covclc
syn keyword FortranFunction dfault
syn keyword FortranFunction dotprd
syn keyword FortranFunction dupdat
syn keyword FortranFunction gqtstp
syn keyword FortranFunction itsmry
syn keyword FortranFunction linvrt
syn keyword FortranFunction litvmu
syn keyword FortranFunction livmul
syn keyword FortranFunction lmstep
syn keyword FortranFunction lsqrt
syn keyword FortranFunction lsvmin
syn keyword FortranFunction ltsqar
syn keyword FortranFunction parchk
syn keyword FortranFunction qapply
syn keyword FortranFunction qrfact
syn keyword FortranFunction reldst
syn keyword FortranFunction rptmul
syn keyword FortranFunction slupdt
syn keyword FortranFunction slvmul
syn keyword FortranFunction stopx
syn keyword FortranFunction vaxpy
syn keyword FortranFunction vcopy
syn keyword FortranFunction vscopy
syn keyword FortranFunction v2norm
let b:num_functions = 32
let b:string = ['imdcon', 'rmdcon', 'errstop', 'nl2sol', 'nl2sno', 'nl2itr', 'assess', 'covclc', 'dfault', 'dotprd', 'dupdat', 'gqtstp', 'itsmry', 'linvrt', 'litvmu', 'livmul', 'lmstep', 'lsqrt', 'lsvmin', 'ltsqar', 'parchk', 'qapply', 'qrfact', 'reldst', 'rptmul', 'slupdt', 'slvmul', 'stopx', 'vaxpy', 'vcopy', 'vscopy', 'v2norm']
let b:type = ['function', 'function', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'function', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'subroutine', 'function', 'subroutine', 'subroutine', 'subroutine', 'function', 'subroutine', 'subroutine', 'subroutine', 'function']
let b:start_line = [9, 26, 55, 99, 655, 834, 1512, 2019, 2402, 2469, 2500, 2537, 3109, 3364, 3402, 3435, 3465, 3910, 3964, 4118, 4154, 4349, 4428, 4638, 4663, 4726, 4765, 4798, 4818, 4831, 4844, 4857]
let b:end_line = [23, 53, 59, 652, 831, 1509, 2016, 2399, 2466, 2497, 2534, 3106, 3361, 3399, 3432, 3462, 3907, 3961, 4115, 4151, 4346, 4425, 4635, 4660, 4723, 4762, 4795, 4815, 4828, 4841, 4854, 4906]
let b:start_line_bracket = [10, 27, 56, 100, 656, 835, 1513, 2020, 2403, 2470, 2501, 2538, 3110, 3365, 3403, 3436, 3466, 3911, 3965, 4119, 4155, 4350, 4429, 4639, 4664, 4727, 4766, 4799, 4819, 4832, 4845, 4858]
syn keyword FortranModule machine_constants
syn keyword FortranModule toms573
syn match FortranModuleVariable display /\%>0l\%<63l\%(%\@<!\<dp\>\)/
syn match FortranModuleVariable display /\%>64l\%<4910l\%(%\@<!\<stop_nl2sol\>\)/
let b:direct_modules_used_list = ['machine_constants', 'machine_constants', 'machine_constants', 'machine_constants']
let b:direct_modules_use_all = ['machine_constants', 'machine_constants', 'machine_constants', 'machine_constants']
let b:direct_modules_only_dict = {}
let b:direct_rename_dict = {}
let b:complete_module_dict = {'interface:imdcon': {'arg:order': ['k'], 'k': ['INTEGER', 'INTENT(in)']}, 'imdcon': ['INTEGER', 'function', 'module:machine_constants', 'nl2sol.f90'], 'dp': ['INTEGER', 'PARAMETER', '=selected_real_kind(precision(1.d0),range(1.d0))', 'variable', 'module:machine_constants', 'nl2sol.f90'], 'interface:errstop': {'routine': ['CHARACTER(*)', 'INTENT(in)'], 'arg:order': ['routine', 'error'], 'error': ['CHARACTER(*)', 'INTENT(in)']}, 'errstop': ['subroutine', 'module:machine_constants', 'nl2sol.f90'], 'machine_constants': ['module', 'nl2sol.f90'], 'interface:rmdcon': {'arg:order': ['k'], 'k': ['INTEGER', 'INTENT(in)']}, 'rmdcon': ['REAL(dp)', 'function', 'module:machine_constants', 'nl2sol.f90']}
syn match FortranFunctionArgument display /\%>8l\%<24l%\@<!\<k\>/
syn match FortranFunctionArgument display /\%>25l\%<54l%\@<!\<k\>/
syn match FortranFunctionArgument display /\%>54l\%<60l%\@<!\<routine\>/
syn match FortranFunctionArgument display /\%>54l\%<60l%\@<!\<error\>/
syn match FortranFunctionArgument display /\%>98l\%<653l%\@<!\<n\>/
syn match FortranFunctionArgument display /\%>98l\%<653l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>98l\%<653l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>98l\%<653l%\@<!\<calcr\>/
syn match FortranFunctionArgument display /\%>98l\%<653l%\@<!\<calcj\>/
syn match FortranFunctionArgument display /\%>98l\%<653l%\@<!\<iv\>/
syn match FortranFunctionArgument display /\%>98l\%<653l%\@<!\<v\>/
syn match FortranFunctionArgument display /\%>98l\%<653l%\@<!\<uiparm\>/
syn match FortranFunctionArgument display /\%>98l\%<653l%\@<!\<urparm\>/
syn match FortranFunctionArgument display /\%>98l\%<653l%\@<!\<ufparm\>/
syn match FortranFunctionArgument display /\%>654l\%<832l%\@<!\<n\>/
syn match FortranFunctionArgument display /\%>654l\%<832l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>654l\%<832l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>654l\%<832l%\@<!\<calcr\>/
syn match FortranFunctionArgument display /\%>654l\%<832l%\@<!\<iv\>/
syn match FortranFunctionArgument display /\%>654l\%<832l%\@<!\<v\>/
syn match FortranFunctionArgument display /\%>654l\%<832l%\@<!\<uiparm\>/
syn match FortranFunctionArgument display /\%>654l\%<832l%\@<!\<urparm\>/
syn match FortranFunctionArgument display /\%>654l\%<832l%\@<!\<ufparm\>/
syn match FortranFunctionArgument display /\%>833l\%<1510l%\@<!\<d\>/
syn match FortranFunctionArgument display /\%>833l\%<1510l%\@<!\<iv\>/
syn match FortranFunctionArgument display /\%>833l\%<1510l%\@<!\<jac\>/
syn match FortranFunctionArgument display /\%>833l\%<1510l%\@<!\<n\>/
syn match FortranFunctionArgument display /\%>833l\%<1510l%\@<!\<nn\>/
syn match FortranFunctionArgument display /\%>833l\%<1510l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>833l\%<1510l%\@<!\<r\>/
syn match FortranFunctionArgument display /\%>833l\%<1510l%\@<!\<v\>/
syn match FortranFunctionArgument display /\%>833l\%<1510l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>1511l\%<2017l%\@<!\<d\>/
syn match FortranFunctionArgument display /\%>1511l\%<2017l%\@<!\<iv\>/
syn match FortranFunctionArgument display /\%>1511l\%<2017l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>1511l\%<2017l%\@<!\<step\>/
syn match FortranFunctionArgument display /\%>1511l\%<2017l%\@<!\<stlstg\>/
syn match FortranFunctionArgument display /\%>1511l\%<2017l%\@<!\<v\>/
syn match FortranFunctionArgument display /\%>1511l\%<2017l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>1511l\%<2017l%\@<!\<x0\>/
syn match FortranFunctionArgument display /\%>2018l\%<2400l%\@<!\<covirc\>/
syn match FortranFunctionArgument display /\%>2018l\%<2400l%\@<!\<d\>/
syn match FortranFunctionArgument display /\%>2018l\%<2400l%\@<!\<iv\>/
syn match FortranFunctionArgument display /\%>2018l\%<2400l%\@<!\<j\>/
syn match FortranFunctionArgument display /\%>2018l\%<2400l%\@<!\<n\>/
syn match FortranFunctionArgument display /\%>2018l\%<2400l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>2018l\%<2400l%\@<!\<r\>/
syn match FortranFunctionArgument display /\%>2018l\%<2400l%\@<!\<v\>/
syn match FortranFunctionArgument display /\%>2018l\%<2400l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>2401l\%<2467l%\@<!\<iv\>/
syn match FortranFunctionArgument display /\%>2401l\%<2467l%\@<!\<v\>/
syn match FortranFunctionArgument display /\%>2468l\%<2498l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>2468l\%<2498l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>2468l\%<2498l%\@<!\<y\>/
syn match FortranFunctionArgument display /\%>2499l\%<2535l%\@<!\<d\>/
syn match FortranFunctionArgument display /\%>2499l\%<2535l%\@<!\<iv\>/
syn match FortranFunctionArgument display /\%>2499l\%<2535l%\@<!\<j\>/
syn match FortranFunctionArgument display /\%>2499l\%<2535l%\@<!\<n\>/
syn match FortranFunctionArgument display /\%>2499l\%<2535l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>2499l\%<2535l%\@<!\<v\>/
syn match FortranFunctionArgument display /\%>2536l\%<3107l%\@<!\<d\>/
syn match FortranFunctionArgument display /\%>2536l\%<3107l%\@<!\<dig\>/
syn match FortranFunctionArgument display /\%>2536l\%<3107l%\@<!\<dihdi\>/
syn match FortranFunctionArgument display /\%>2536l\%<3107l%\@<!\<ka\>/
syn match FortranFunctionArgument display /\%>2536l\%<3107l%\@<!\<l\>/
syn match FortranFunctionArgument display /\%>2536l\%<3107l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>2536l\%<3107l%\@<!\<step\>/
syn match FortranFunctionArgument display /\%>2536l\%<3107l%\@<!\<v\>/
syn match FortranFunctionArgument display /\%>2536l\%<3107l%\@<!\<w\>/
syn match FortranFunctionArgument display /\%>3108l\%<3362l%\@<!\<d\>/
syn match FortranFunctionArgument display /\%>3108l\%<3362l%\@<!\<iv\>/
syn match FortranFunctionArgument display /\%>3108l\%<3362l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>3108l\%<3362l%\@<!\<v\>/
syn match FortranFunctionArgument display /\%>3108l\%<3362l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>3363l\%<3400l%\@<!\<n\>/
syn match FortranFunctionArgument display /\%>3363l\%<3400l%\@<!\<lin\>/
syn match FortranFunctionArgument display /\%>3363l\%<3400l%\@<!\<l\>/
syn match FortranFunctionArgument display /\%>3401l\%<3433l%\@<!\<n\>/
syn match FortranFunctionArgument display /\%>3401l\%<3433l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>3401l\%<3433l%\@<!\<l\>/
syn match FortranFunctionArgument display /\%>3401l\%<3433l%\@<!\<y\>/
syn match FortranFunctionArgument display /\%>3434l\%<3463l%\@<!\<n\>/
syn match FortranFunctionArgument display /\%>3434l\%<3463l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>3434l\%<3463l%\@<!\<l\>/
syn match FortranFunctionArgument display /\%>3434l\%<3463l%\@<!\<y\>/
syn match FortranFunctionArgument display /\%>3464l\%<3908l%\@<!\<d\>/
syn match FortranFunctionArgument display /\%>3464l\%<3908l%\@<!\<g\>/
syn match FortranFunctionArgument display /\%>3464l\%<3908l%\@<!\<ierr\>/
syn match FortranFunctionArgument display /\%>3464l\%<3908l%\@<!\<ipivot\>/
syn match FortranFunctionArgument display /\%>3464l\%<3908l%\@<!\<ka\>/
syn match FortranFunctionArgument display /\%>3464l\%<3908l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>3464l\%<3908l%\@<!\<qtr\>/
syn match FortranFunctionArgument display /\%>3464l\%<3908l%\@<!\<r\>/
syn match FortranFunctionArgument display /\%>3464l\%<3908l%\@<!\<step\>/
syn match FortranFunctionArgument display /\%>3464l\%<3908l%\@<!\<v\>/
syn match FortranFunctionArgument display /\%>3464l\%<3908l%\@<!\<w\>/
syn match FortranFunctionArgument display /\%>3909l\%<3962l%\@<!\<n1\>/
syn match FortranFunctionArgument display /\%>3909l\%<3962l%\@<!\<n\>/
syn match FortranFunctionArgument display /\%>3909l\%<3962l%\@<!\<l\>/
syn match FortranFunctionArgument display /\%>3909l\%<3962l%\@<!\<a\>/
syn match FortranFunctionArgument display /\%>3909l\%<3962l%\@<!\<irc\>/
syn match FortranFunctionArgument display /\%>3963l\%<4116l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>3963l\%<4116l%\@<!\<l\>/
syn match FortranFunctionArgument display /\%>3963l\%<4116l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>3963l\%<4116l%\@<!\<y\>/
syn match FortranFunctionArgument display /\%>3963l\%<4116l%\@<!\<fn_val\>/
syn match FortranFunctionArgument display /\%>4117l\%<4152l%\@<!\<n\>/
syn match FortranFunctionArgument display /\%>4117l\%<4152l%\@<!\<a\>/
syn match FortranFunctionArgument display /\%>4117l\%<4152l%\@<!\<l\>/
syn match FortranFunctionArgument display /\%>4153l\%<4347l%\@<!\<iv\>/
syn match FortranFunctionArgument display /\%>4153l\%<4347l%\@<!\<n\>/
syn match FortranFunctionArgument display /\%>4153l\%<4347l%\@<!\<nn\>/
syn match FortranFunctionArgument display /\%>4153l\%<4347l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>4153l\%<4347l%\@<!\<v\>/
syn match FortranFunctionArgument display /\%>4348l\%<4426l%\@<!\<n\>/
syn match FortranFunctionArgument display /\%>4348l\%<4426l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>4348l\%<4426l%\@<!\<j\>/
syn match FortranFunctionArgument display /\%>4348l\%<4426l%\@<!\<r\>/
syn match FortranFunctionArgument display /\%>4348l\%<4426l%\@<!\<ierr\>/
syn match FortranFunctionArgument display /\%>4427l\%<4636l%\@<!\<m\>/
syn match FortranFunctionArgument display /\%>4427l\%<4636l%\@<!\<n\>/
syn match FortranFunctionArgument display /\%>4427l\%<4636l%\@<!\<qr\>/
syn match FortranFunctionArgument display /\%>4427l\%<4636l%\@<!\<alpha\>/
syn match FortranFunctionArgument display /\%>4427l\%<4636l%\@<!\<ipivot\>/
syn match FortranFunctionArgument display /\%>4427l\%<4636l%\@<!\<ierr\>/
syn match FortranFunctionArgument display /\%>4427l\%<4636l%\@<!\<nopivk\>/
syn match FortranFunctionArgument display /\%>4427l\%<4636l%\@<!\<sum\>/
syn match FortranFunctionArgument display /\%>4637l\%<4661l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>4637l\%<4661l%\@<!\<d\>/
syn match FortranFunctionArgument display /\%>4637l\%<4661l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>4637l\%<4661l%\@<!\<x0\>/
syn match FortranFunctionArgument display /\%>4662l\%<4724l%\@<!\<func\>/
syn match FortranFunctionArgument display /\%>4662l\%<4724l%\@<!\<ipivot\>/
syn match FortranFunctionArgument display /\%>4662l\%<4724l%\@<!\<j\>/
syn match FortranFunctionArgument display /\%>4662l\%<4724l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>4662l\%<4724l%\@<!\<rd\>/
syn match FortranFunctionArgument display /\%>4662l\%<4724l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>4662l\%<4724l%\@<!\<y\>/
syn match FortranFunctionArgument display /\%>4662l\%<4724l%\@<!\<z\>/
syn match FortranFunctionArgument display /\%>4725l\%<4763l%\@<!\<a\>/
syn match FortranFunctionArgument display /\%>4725l\%<4763l%\@<!\<cosmin\>/
syn match FortranFunctionArgument display /\%>4725l\%<4763l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>4725l\%<4763l%\@<!\<size\>/
syn match FortranFunctionArgument display /\%>4725l\%<4763l%\@<!\<step\>/
syn match FortranFunctionArgument display /\%>4725l\%<4763l%\@<!\<u\>/
syn match FortranFunctionArgument display /\%>4725l\%<4763l%\@<!\<w\>/
syn match FortranFunctionArgument display /\%>4725l\%<4763l%\@<!\<wchmtd\>/
syn match FortranFunctionArgument display /\%>4725l\%<4763l%\@<!\<wscale\>/
syn match FortranFunctionArgument display /\%>4725l\%<4763l%\@<!\<y\>/
syn match FortranFunctionArgument display /\%>4764l\%<4796l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>4764l\%<4796l%\@<!\<y\>/
syn match FortranFunctionArgument display /\%>4764l\%<4796l%\@<!\<s\>/
syn match FortranFunctionArgument display /\%>4764l\%<4796l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>4817l\%<4829l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>4817l\%<4829l%\@<!\<w\>/
syn match FortranFunctionArgument display /\%>4817l\%<4829l%\@<!\<a\>/
syn match FortranFunctionArgument display /\%>4817l\%<4829l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>4817l\%<4829l%\@<!\<y\>/
syn match FortranFunctionArgument display /\%>4830l\%<4842l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>4830l\%<4842l%\@<!\<y\>/
syn match FortranFunctionArgument display /\%>4830l\%<4842l%\@<!\<x\>/
syn match FortranFunctionArgument display /\%>4843l\%<4855l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>4843l\%<4855l%\@<!\<y\>/
syn match FortranFunctionArgument display /\%>4843l\%<4855l%\@<!\<s\>/
syn match FortranFunctionArgument display /\%>4856l\%<4907l%\@<!\<p\>/
syn match FortranFunctionArgument display /\%>4856l\%<4907l%\@<!\<x\>/
syn match FortranUnusedVariable display /\%>54l\%<60l%\@<!\<routine\>/
syn match FortranUnusedVariable display /\%>54l\%<60l%\@<!\<error\>/
syn match FortranUnusedVariable display /\%>4153l\%<4347l%\@<!\<vn\>/