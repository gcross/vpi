<?xml version="1.0" encoding="utf-8"?>
<?xml-stylesheet ekr_test?>
<leo_file>
<leo_header file_format="2" tnodes="0" max_tnode_index="0" clone_windows="0"/>
<globals body_outline_ratio="0.5">
	<global_window_position top="22" left="10" height="852" width="804"/>
	<global_log_window_position top="0" left="0" height="0" width="0"/>
</globals>
<preferences/>
<find_panel_settings/>
<vnodes>
<v t="gcross.20090623152316.1" str_leo_pos="1,3,4,2"><vh>@chapters</vh></v>
<v t="gcross.20090623152316.2" a="E"><vh>@path src</vh>
<v t="gcross.20090623152316.4"><vh>&lt;&lt; Choice of trial / potential functions &gt;&gt;</vh></v>
<v t="gcross.20090623152316.3"><vh>@thin vpi_main.f90</vh></v>
<v t="gcross.20090623152316.34" a="E"
expanded="gcross.20090623152316.56,gcross.20090623184136.3812,"><vh>@thin vpi_potential.f90</vh></v>
<v t="gcross.20090623152316.73" a="E"
expanded="gcross.20090623152316.80,"><vh>@thin vpi_sample.f90</vh></v>
<v t="gcross.20090623152316.91"><vh>@thin vpi_trial_func.f90</vh></v>
</v>
</vnodes>
<tnodes>
<t tx="gcross.20090623152316.1"></t>
<t tx="gcross.20090623152316.2"></t>
<t tx="gcross.20090623152316.4">use vpi_potential, Usp_func =&gt; vpi_Usp_HOz, &amp;
                   Uij_func =&gt; vpi_Uij_NULL, &amp;
                   gUsp_func =&gt; vpi_gUsp_HOz, &amp;
                   gUij_func =&gt; vpi_gUij_NULL, &amp;
                   Uij_rot_func =&gt; vpi_RUij_undefined
use vpi_trial_func, tfunc =&gt; hoz_tfun, &amp;
                    jas_tfun =&gt; independent_tfun, &amp;
                    grad_lap_sp_tfun =&gt; numeric_grad_lap_spf, &amp;
                    grad_lap_jas_tfun =&gt; grad_lap_independent_tfun
</t>
</tnodes>
</leo_file>
