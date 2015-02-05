This code has been adopted from the original WHAMP code version
to allow model input from files and a few more output parameters.
The original WHAMP version  can be downloaded from
http://www.tp.umu.se/forskning/space/WHAMP/
It is also included under WHAMP_code directory with
the name "WHAMP.tar.gz_original"

The output parameters from WHAMP (can be specified
running WHAMP) are

============ output ===========
        e       (ex,ey,ez)
        b       (bx,by,bz)
        f       frequency <real,imaginery>
        g       group velocity
        h       |e|/|b| [mV/nT]
        l       |bp|/|bz|
        m       Im[bx]/Re[by]
        n       ellipticity bx/by
        o       ellipticity general
        p       k perpendicular
        s       (sp, sz)   spatial growth
        u       energy ration between the total wave energy and energy in electric field
        v       Poynting flux (in uW/m2 for <E^2>=0.5(mV/m)^2)
        z       k paralel
        x       phase of bz against bx (/+/ means bz is in front of bx)
        y       energy density and flux of each component
================================

the input model file looks like

============ Model files ===========
n(1)  n(2)  ...  n(10)   /per m3/
t(1)  t(2)  ...  t(10)   / keV, T_par /
d(1)  d(2)  ...  d(10)   / loss cone parameter, default 1.0 (no loss cone) /
a(1)  a(2)  ...  a(10)   / t_perp/t_par, default 1.0 /
b(1)  b(2)  ...  b(10)   / default 0, i.e. no loss cone/
ass(1)  ass(2)  ...  ass(10)    / 0-electrons, 1-protons, 16-oxygen /
vd(1)  vd(2)  ...  vd(10)   / v_drift/v_term /
fce / electron gyrofrequency in kHz /
pzl     / 1 - log scale, 0 - linear scale /
=====================================
Plasma with B=100nT, n=1cm-3, there are only
oxygen ions (Tperp=50eV, Tpar=10eV, vd=1*vt_par)
and 1eV electrons.
============ Example_Input_File would have 9 lines ===========
1.e6 1.e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
.001 .001 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
1.   1.   1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
5.   1.   1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
0.0  0.0  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
16.  0.0  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
1.   0.0  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
2.8
0


Making MEX file in MATLAB
=========================
>mex libwhamp.a mexwhamp.F

