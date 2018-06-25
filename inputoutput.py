#!/usr/bin/env python


input_vars = [
    # Type, Name, Dim, , Length
    ['character(16)', 'sysname', '', '"untitled"', 1],
    ['integer', 'nx1_m', '', -500, 1],
    ['integer', 'ny1_m', '', 1, 1],
    ['integer', 'nz1_m', '', 1, 1],
    ['integer', 'nx2_m', '', +500, 1],
    ['integer', 'ny2_m', '', 1, 1],
    ['integer', 'nz2_m', '', 1, 1],
    ['integer', 'nt', '', 1000, 1],
    ['real(8)', 'hx_m', '', '250d0', 1],
    ['real(8)', 'hy_m', '', '250d0', 1],
    ['real(8)', 'hz_m', '', '250d0', 1],
    ['real(8)', 'dt', '', 1, 1],
    ['real(8)', 'ac_0', '', '1d0', 1],
    ['real(8)', 'epdir', '(1:3)', '(/0d0,0d0,1d0/)', 3],
    ['real(8)', 't_pulse', '', '440d0', 1],
    ['real(8)', 'omega', '', '0.057d0', 1],
    ['real(8)', 'omega_l', '', '1d1', 1],
    ['real(8)', 'gamma_l', '', '0d0', 1],
    ['real(8)', 'chi_l0', '', '1d0', 1],
    ['logical', "out_ac_bin", '', '.false.', 1],
    ['logical', "out_ac_out", '', '.true.', 1],
    ['integer', "ac_bin_step", '', '10', 1],
    ['integer', "ac_out_step", '', '1000', 1],
    ['real(8)', "angle", '', '0', 1],
    ['logical', "inp_bc", '', '.true.', 1],
]
























template = r"""
module inputoutput
  implicit none

{DEF_VARS}

contains



  subroutine read_input()
    implicit none

{DEF_NML}

{INIT_VARS}

    read (*, nml=input)

    return
  end subroutine read_input



  subroutine var_dump()
    implicit none
    
    print '(a)', "################ VAR_DUMP ################"
{VAR_DUMP}
    print '(a)', "##########################################"
    return
  end subroutine var_dump



end module inputoutput
"""

def_vars_list = []
def_nml_list = []
init_vars_list = []
var_dump_list = []

vfmt = {
    'integer': 'i6',
    'real(8)': 'es23.15e3',
    'character(16)': 'a',
    'logical': 'l1',
}

def add_indent(code, level=0):
    return '\n'.join([
        (' '*level + line)
        for line in code.splitlines()
    ])

for vtype, vname, vdim, vdefval, vlen in input_vars:
    def_vars_list += ['{VTYPE} :: {VNAME}{VDIM}'.format(
        VTYPE=vtype, VNAME=vname, VDIM=vdim
    )]
    def_nml_list += [vname]
    init_vars_list += ['{VNAME}={VDEFVAL}'.format(
        VNAME=vname, VDEFVAL=vdefval
    )]
    var_dump_list += ["""print '("# {VNAME}", {VLEN}(1X, {VFMT}))', {VNAME}""".format(
        VNAME=vname, VLEN=vlen, VFMT=vfmt[vtype]
    )
    ]

with open(__file__.replace('.py', '') + '.f90', 'w') as fh:
    fh.write(template.format(
        DEF_VARS = add_indent("\n".join(def_vars_list),2),
        DEF_NML = add_indent("namelist/input/ &\n& "  + ", &\n& ".join(def_nml_list),4),
        INIT_VARS = add_indent("\n".join(init_vars_list),4),
        VAR_DUMP = add_indent("\n".join(var_dump_list),4),
    ))
    print("Wrote %s" % fh.name)
