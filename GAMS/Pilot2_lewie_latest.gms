$title PILOT: FIRST GENERIC LEWIE MODEL

* This model is a generic LEWIE.  The data is randomely drawn from known distributions of variables.
* The model is set up to solve repeatedly using consecutive random draws of input variables and parameters.
* The multiple runs of the model allow to compute mean/stdev for the outcome variables.


*----------------------------------------------------------------------------------------------
* choose the number of draws (the second number)
set draw /dr1*dr3/ ;
*----------------------------------------------------------------------------------------------

* A few useful gams options
option limrow=10 ;
option limcol=10 ;
*$onsymlist
$onsymxref
* unstar the following only if you don't have a PATH licence
*option mcp = miles;


* ----------------------------------------------------------------------------------------------
* ------------------------------- GENERIC MODEL SPECIFICATION -----------------------------------
* ----------------------------------------------------------------------------------------------

* First define all sets
* -----------------------
sets
* goods
     g         commodities /g1*g3/

* factors
     f         inputs /f1*f4/
     fk(f)     fixed factors /f1,f2/
     ft(f)     tradable inputs /f3,f4/

* households
     h         households /h1*h4/

* exogenous accounts **** not in this version yet!
*     n         exogenous /zoi, row, gov/
;

alias (g,gg)
      (h,hh)
      (f,fa) ;


* Now variables and parameters
* ---------------------------------
nonnegative variables
* production
     QP(g,h)   quantity produced of a good by a household
     FD(g,f,h) factor demand of f in production of g
*     ID(g,gg,h) intermediate demand for production of g  **** not in this version yet!
     QFS(f)   quantity of factor supply (for traded factors)

* consumption
     QC(g,h)   quantity consumed of g by h
     Y(h)      household income

* values
     PM(g)     market price of a good
     WK(g,f,h) value of fixed factors
     WT(f)     value of tradable factors
;

variables
* trade
     QX(g)     net quantity exported of g
     USELESSVAR   fake variable for the nlp versions of the model
;
USELESSVAR.l = 1 ;

parameters
*Production -
*CES
     aces(g,h) shift parameter in the CES production function
     shces(g,f,h) factor share parameters in the CES
     rhoces(g,h)  substitution elasticity in CES production functions
     sigces(g,h)  sigma is 1 over (1 minus rho)

*Cobb-douglas
     acobb(g,h) production shift parameter for the CD
     shcobb(g,f,h) factor share parameter for the CD

*Consumption
     alpha(g,h) consumption share parameters in the LES
     cmin(g,h)  minimal consumption in the LES
     exinc(h)   exogenous income of household

* factor endowments for fixed factors
     fixfac(g,f,h) fixed factors

* Factor supply
     se(f) factor supply elasticity for tradable factors
;


* Initial values for parameters
parameter
     qp0(g,h)   initial quantity produced
     qc0(g,h)   initial quantity consumed
     fd0(g,f,h) initial factor demands
     qfs0(f)    initial tradable factor supply
     wk0(g,f,h) initial factor values
     wt0(f)     initial factor values
     qx0(g)     initial net exports
     y0(h)      initial household income
     pm0(g)     initial market price
;

Equations
* production
     EQ_FDCOBB(g,f,h)    factor demands cobb douglas
     EQ_FDLEON(g,f,h)    factor demands CES
     EQ_FDCES(g,f,h)     factor demands CES other version
     EQ_FDCES2(g,f,fa,h) factor demands CES other version
     EQ_QPCOBB(g,h)      quantity produced cobb douglas
     EQ_QPLEON(g,h)      quantity produced leontieff
     EQ_QPCES(g,h)       quantity produced CES
     EQ_QPCES2(g,h)      quantity produced CES2
     EQ_ZEROPROF(g,h)    zero profit condition
     EQ_FACSUP(f)        factor supply for tradable factors

* consumption
     EQ_QC(g,h)          quantity consumed

* income
     EQ_Y(h)             full income constraint for the household
* market clearing
     EQ_MKT(g)           market clearing in the economy
     EQ_FMKT1(g,f,h)     fixed factors clearing
     EQ_FMKT2(f)         tradable factors clearing
;


* PRODUCTION BLOCK
* Cobb douglas form
EQ_QPCOBB(g,h)..   QP(g,h) =E= acobb(g,h)*prod(f,FD(g,f,h)**(shcobb(g,f,h)))
;
EQ_FDCOBB(g,f,h)..
     (WK(g,f,h)*FD(g,f,h)- PM(g)*QP(g,h)*shcobb(g,f,h))$fk(f)
    +(WT(f)*FD(g,f,h)- PM(g)*QP(g,h)*shcobb(g,f,h))$ft(f)
    =E= 0
;

* CES form
EQ_QPCES(g,h)..
     QP(g,h) =E= qp0(g,h)*(sum(f,shces(g,f,h)*(FD(g,f,h)/fd0(g,f,h))**(rhoces(g,h))))**(1/rhoces(g,h))
;

EQ_FDCES(g,f,h)..
       [FD(g,f,h)/fd0(g,f,h) - (QP(g,h)/qp0(g,h))*(PM(g)/1*wk0(g,f,h)/WK(g,f,h))**(sigces(g,h))]$fk(f)
      +[FD(g,f,h)/fd0(g,f,h) - (QP(g,h)/qp0(g,h))*(PM(g)/1*wt0(f)/WT(f))**(sigces(g,h))]$ft(f)
      =E= 0
;

* CONSUMPTION AND INCOME
EQ_QC(g,h)..
     QC(g,h) =E= alpha(g,h)/PM(g)*[Y(h)-sum(gg, PM(gg)*cmin(gg,h))] + cmin(g,h)
;

* Full income (value of factor endowments)
EQ_Y(h)..
     Y(h) =E= sum((g,fk),WK(g,fk,h)*FD(g,fk,h)) + exinc(h)
;

* FACTOR MARKETS
* market clearing
EQ_MKT(g)..
     sum(h,QP(g,h)-QC(g,h)) =E= QX(g)
;

* market clearing for factors
EQ_FMKT1(g,f,h)$fk(f)..
        FD(g,f,h) =E= fixfac(g,f,h)
;

* fixed factor market clearing (=rent determination)
EQ_FMKT2(ft)..
        sum((g,h),FD(g,ft,h)) =E= QFS(ft)
;

* Factor supply elasticity
EQ_FACSUP(ft)..
     QFS(ft)/qfs0(ft) =e= (WT(ft)/wt0(ft))**se(ft) ;
;

*-------------------------------------------------------------------------------------------------
*--------------------------------------- ALTERNATIVE MODELS --------------------------------------
*-------------------------------------------------------------------------------------------------

model genCD Model with Cobb Douglas production /
      EQ_QPCOBB.QP
      EQ_QC.QC
      EQ_FDCOBB.FD
      EQ_MKT.QX
      EQ_Y.Y
      EQ_FMKT1.WK
      EQ_FMKT2.QFS
      EQ_FACSUP.WT
/;

model genCES Model with CES production /
      EQ_QPCES.QP
      EQ_QC.QC
      EQ_FDCES.FD
      EQ_MKT.QX
      EQ_Y.Y
      EQ_FMKT1.WK
      EQ_FMKT2.QFS
      EQ_FACSUP.WT
/;


* ----------------------------------------------------------------------------------------------
* ------------------------------- END OF MODEL SPECIFICATION -----------------------------------
* ----------------------------------------------------------------------------------------------

* ----------------------------------------------------------------------------------------------
* ------------------------------- START OF DATA CALIBRATION  -----------------------------------
* ----------------------------------------------------------------------------------------------

* We don't draw each iteration of the dataset while looping. Instead, we create a parameter with multiple draws, which
* can then be read as input values for each loop.
* Thos parameters will all have a "_dr" suffix, and are indexed by the "draw" set.


* 1. parameters for which we actually draw a value
*------------------------------------------------
* NOTE: ONLY 4 PARAMETERS ARE TRULY RANDOM - OTHERS ARE UNIQUELY DETERMINED FROM THOSE 4
parameter
     sigces_dr(g,h,draw) draw from sigces distribution
     fd_dr(g,f,h,draw) draw from the fd distribution
     alpha_dr(g,h,draw) consumption alphas in each draw
     se_dr(f,draw)  supply elasticity of tradable factors

;
* choose you randomisation seed if you want repetible and comparable outputs
* choose the execseed line instead if you want true randomness
option seed = 3000;
*execseed = gmillisec(jnow);

* elasticity of factor substitution in CES in each draw
sigces_dr(g,h,draw) = uniform(0.3,0.5) ;

* Factor demand in CES in each draw
* this is very generic:
fd_dr(g,f,h,draw) = uniform(0,0.5);
* But we could have different distributions for different goods/households:
* fd_dr(g,f,"h2",draw) = uniform(0.5,1);
* fd_dr("g1",f,"h3",draw) = uniform(0.3,1.5) ;

* consumption shares in each draw
parameter temp(g,h,draw) ;
temp(g,h,draw) = uniform(0,1);    display temp ;
alpha_dr(g,h,draw) = temp(g,h,draw) / sum(gg,temp(gg,h,draw));

* supply elasticity in each draw
se_dr(ft,draw) = uniform(0.2,0.8);

display sigces_dr, fd_dr, alpha_dr, se_dr;


* 2. parameters which we just compute from values drawn randomely
* -------------------------------------------------------------
parameter
     rhoces_dr(g,h,draw) CES rho comping from the draw of sigma
     qp_dr(g,h,draw) QP in each draw
     shcobb_dr(g,f,h,draw) cobb douglas share param in each draw
     acobb_dr(g,h,draw) cobb douglas shift parameter in each draw
     shces_dr(g,f,h,draw) CES share parameter in each draw
     aces_dr(g,h,draw) CES shift parameter in each draw
     shces2_dr(g,f,h,draw) CES share parameter in each draw
     aces2_dr(g,h,draw) CES shift parameter in each draw
     fixfac_dr(g,f,h,draw) fixed factors in each draw
     qfs_dr(f,draw) quantity of tradable factors supplied in each draw
     y_dr(h,draw) household income in each draw
     exinc_dr(h,draw) exogenous income in each draw
     qc_dr(g,h,draw) quantities consumed in each draw
     qx_dr(g,draw) quantitites exported in each draw
;

* rho supposed to be in ]-inf,0[U]0,1[
* sigma supposed to be in ]0,1[U]1,+inf[
* sigma is more easily interpreted (factor elasticity of substitution)
rhoces_dr(g,h,draw) = (sigces_dr(g,h,draw)-1)/sigces_dr(g,h,draw);
display sigces_dr, rhoces_dr;

* factor prices all =1 in the baseline
wk0(g,fk,h)=1;
wt0(ft)=1;
display wk0, wt0 ;

* in the baseline, factor prices are one, so qp is the sum of fd:
qp_dr(g,h,draw)=sum(f,fd_dr(g,f,h,draw));
display qp_dr ;

* cobb douglas parameters (determined from the drawn factor demands):
shcobb_dr(g,f,h,draw) = fd_dr(g,f,h,draw)/sum(fa,fd_dr(g,fa,h,draw));
acobb_dr(g,h,draw) = qp_dr(g,h,draw)/prod(fa,fd_dr(g,fa,h,draw)**shcobb_dr(g,fa,h,draw)) ;
display shcobb_dr, acobb_dr ;

* CES parameters (determined from the drawn factor demands):
* this calibration method is from Tom Rutherfords lecture notes on CES, in "calibrated share form"
* http://www.gamsworld.eu/mpsge/debreu/ces.pdf
* share parameters (assuming factor prices are all initially 1, otherwise it changes)
shces_dr(g,f,h,draw) = fd_dr(g,f,h,draw)/ (sum(fa,fd_dr(g,fa,h,draw))) ;
* NOTE: there is no shift parameter in the calibrated share form
display fd_dr, shces_dr;

fixfac_dr(g,fk,h,draw) = fd_dr(g,fk,h,draw) ;
qfs_dr(ft,draw) = sum((g,h),fd_dr(g,ft,h,draw)) ;
display fixfac_dr, qfs_dr ;

* income and consumption
* for now exogenous income at 0
exinc_dr(h,draw) = 0 ;
* initial commodity prices all =1
pm0(g) = 1;
* income is determined from value of household assets (with rent/wage = 1)
y_dr(h,draw) = sum((g,fk),wk0(g,fk,h)*fd_dr(g,fk,h,draw)) + exinc_dr(h,draw);
display y_dr ;

* quantity consumed is directly determined from y and consumption shares
* (no need to draw qc, just the shares)
qc_dr(g,h,draw)=y_dr(h,draw)*alpha_dr(g,h,draw)/pm0(g) ;
display qc_dr;
* for now let's keep cmin at zero
cmin(g,h)=0 ;

* initialise qp0 with zero profit condition under price=1:
qp_dr(g,h,draw) = sum(f,fd_dr(g,f,h,draw)*(wk0(g,f,h)$fk(f) + wt0(f)$ft(f))) ;
display qc_dr, qp_dr ;

* net exports
qx_dr(g,draw) = sum(h,qp_dr(g,h,draw)-qc_dr(g,h,draw)) ;
display qx_dr;

* ----------------------------------------------------------------------------------------------
* --------------------------------- END OF DATA CALIBRATION  -----------------------------------
* ----------------------------------------------------------------------------------------------



*---------------------------------------------------------------------------------------------------
*---------------------------------------------------------------------------------------------------
*--------------------------------- NOW THE LOOP ----------------------------------------------------
*---------------------------------------------------------------------------------------------------
*---------------------------------------------------------------------------------------------------

set mod models /cobb, ces/ ;

* parameters to store output in each loop:
parameter
fd1(g,f,h,draw,mod) calibrated fd
qfs1(f,draw,mod)    calibrated qfs
qp1(g,h,draw,mod)   calibrated qp
qc1(g,h,draw,mod)   calibrated qc
qx1(g,draw,mod)     calibrated qx
pm1(g,draw,mod)     calibrated pm
wt1(f,draw,mod)     calibrated wt
wk1(g,f,h,draw,mod) calibrated wk
y1(h,draw,mod)      calibrated y

fd2(g,f,h,draw,mod) simulated fd
qfs2(f,draw,mod)    simulated qfs
qp2(g,h,draw,mod)   simulated qp
qc2(g,h,draw,mod)   simulated qc
qx2(g,draw,mod)     simulated qx
pm2(g,draw,mod)     simulated pm
wt2(f,draw,mod)     simulated wt
wk2(g,f,h,draw,mod) simulated wk
y2(h,draw,mod)      simulated y

fdD(g,f,h,draw,mod) delta fd
qfsD(f,draw,mod)    delta qfs
qpD(g,h,draw,mod)   delta qp
qcD(g,h,draw,mod)   delta qc
qxD(g,draw,mod)     delta qx
pmD(g,draw,mod)     delta pm
wtD(f,draw,mod)     delta wt
wkD(g,f,h,draw,mod) delta wk
yD(h,draw,mod)      delta y

fdPC(g,f,h,draw,mod) percent change fd
qfsPC(f,draw,mod)    percent change qfs
qpPC(g,h,draw,mod)   percent change qp
qcPC(g,h,draw,mod)   percent change qc
qxPC(g,draw,mod)     percent change qx
pmPC(g,draw,mod)     percent change pm
wtPC(f,draw,mod)     percent change wt
wkPC(g,f,h,draw,mod) percent change wk
yPC(h,draw,mod)      percent change y
;



loop(draw,
* Read the parameter values for this iteration of the loop
     acobb(g,h) = acobb_dr(g,h,draw) ;
     shcobb(g,f,h) = shcobb_dr(g,f,h,draw) ;
     shces(g,f,h) = shces_dr(g,f,h,draw);
     rhoces(g,h) = rhoces_dr(g,h,draw);
     sigces(g,h) = sigces_dr(g,h,draw);
     alpha(g,h) = alpha_dr(g,h,draw);
     exinc(h)   = exinc_dr(h,draw);
     fixfac(g,fk,h) =fixfac_dr(g,fk,h,draw);
     se(ft) = se_dr(ft,draw);

* Read initial values for this run of the model:
     FD.l(g,f,h) = fd_dr(g,f,h,draw) ;
     QP.l(g,h) = qp_dr(g,h,draw);
     QC.l(g,h) = qc_dr(g,h,draw);
     QX.l(g) = qx_dr(g,draw);
     WK.l(g,fk,h) = 1 ;
     WT.l(ft) = 1 ;
     PM.l(g) = 1 ;
     QFS.l(ft) = qfs_dr(ft,draw) ;
     Y.l(h) = y_dr(h,draw);

* initialise base values for CES functions in calibrated share form
     qp0(g,h) = qp_dr(g,h,draw) ;
     fd0(g,f,h) = fd_dr(g,f,h,draw) ;
     qfs0(f)   = qfs_dr(f,draw);
*    wk0 and wt0 are already set to 1

* tradable input price determined by elasticity - no need to fix them anymore
     PM.fx(g) =1 ;
     display FD.l, QP.l, QC.l, QX.l, PM.l, WT.l, WK.l, Y.l ;

*---------------------------------
* CALIBRATIONS
*---------------------------------

* 1) Cobb douglas
*-----------------
     solve genCD using mcp ;
     display genCD.modelstat ;
     abort$(not genCD.modelstat eq 1) "Cobb-douglas model calibration failed"
     display "pre-CD" ;
     display FD.l, QP.l, QC.l, QX.l, PM.l, WT.l, WK.l, Y.l ;
     fd1(g,f,h,draw,"cobb") = FD.l(g,f,h);
     qp1(g,h,draw,"cobb") = QP.l(g,h);
     qc1(g,h,draw,"cobb") = QC.l(g,h);
     qx1(g,draw,"cobb") = QX.l(g);
     pm1(g,draw,"cobb") = PM.l(g);
     wt1(f,draw,"cobb") = WT.l(f);
     wk1(g,f,h,draw,"cobb") =WK.l(g,f,h);
     qfs1(f,draw,"cobb") =QFS.l(f);
     y1(h,draw,"cobb") = Y.l(h);

     display fd1, qp1, qc1, qx1, pm1, wt1, wk1, qfs1, y1 ;

* 2) CES
*---------

* reset initial values
     FD.l(g,f,h) = fd_dr(g,f,h,draw) ;
     QP.l(g,h) = qp_dr(g,h,draw);
     QC.l(g,h) = qc_dr(g,h,draw);
     QX.l(g) = qx_dr(g,draw);
     WK.l(g,fk,h) = 1 ;
     WT.l(ft) = 1 ;
     PM.l(g) = 1 ;
     QFS.l(ft) = qfs_dr(ft,draw) ;
     Y.l(h) = y_dr(h,draw);

     solve genCES using mcp ;
     display genCES.modelstat ;
     abort$(not genCES.modelstat eq 1) "CES model calibration failed"
     display "pre-CES"
     display FD.l, QP.l, QC.l, QX.l, PM.l, WT.l, WK.l, Y.l ;

     fd1(g,f,h,draw,"ces") = FD.l(g,f,h);
     qp1(g,h,draw,"ces") = QP.l(g,h);
     qc1(g,h,draw,"ces") = QC.l(g,h);
     qx1(g,draw,"ces") = QX.l(g);
     pm1(g,draw,"ces") = PM.l(g);
     wt1(f,draw,"ces") = WT.l(f);
     wk1(g,f,h,draw,"ces") =WK.l(g,f,h);
     qfs1(f,draw,"ces") =QFS.l(f);
     y1(h,draw,"ces") = Y.l(h);

     display fd1, qp1, qc1, qx1, pm1, wt1, wk1, qfs1 ;

* This is a useful check: check that the calibrations are the same
*$ontext
     loop((g,f,h,mod),
     ABORT$(abs(fd1(g,f,h,draw,mod)-fd1(g,f,h,draw,"cobb")) gt 0.000001) "fd not equal" ;
     ABORT$(abs(qp1(g,h,draw,mod)-qp1(g,h,draw,"cobb")) gt 0.000001) "qp not equal";
     ABORT$(abs(qc1(g,h,draw,mod)-qc1(g,h,draw,"cobb")) gt 0.000001) "qc not equal" ;
     ABORT$(abs(qx1(g,draw,mod)-qx1(g,draw,"cobb")) gt 0.000001)  "qx not equal" ;
     ABORT$(abs(wt1(f,draw,mod)-wt1(f,draw,"cobb")) gt 0.000001) "wt not equal" ;
     ABORT$(abs(wk1(g,f,h,draw,mod)-wk1(g,f,h,draw,"cobb"))gt 0.000001) "wk not equal" ;
     ABORT$(abs(qfs1(f,draw,mod)-qfs1(f,draw,"cobb")) gt 0.000001) "qfs not equal" ;
     );
*$offtext

*--------------------------------------
* SIMULATIONS
*--------------------------------------

* shock the system
*--------------------
     PM.fx("g1") = 1.1  ;

* 1) solve again with Cobb Douglas
*----------------------------------
     solve genCD using mcp ;
     abort$(genCD.modelstat ne 1) "cobb douglas simulation fails" ;
     display "post CD"  ;
     display FD.l, QP.l, QC.l, QX.l, PM.l, WT.l, WK.l, Y.l ;

* Output results
*-----------------
     fd2(g,f,h,draw,"cobb") = FD.l(g,f,h);
     qp2(g,h,draw,"cobb") = QP.l(g,h);
     qc2(g,h,draw,"cobb") = QC.l(g,h);
     qx2(g,draw,"cobb") = QX.l(g);
     pm2(g,draw,"cobb") = PM.l(g);
     wt2(f,draw,"cobb") = WT.l(f);
     wk2(g,f,h,draw,"cobb") =WK.l(g,f,h);
     qfs2(f,draw,"cobb") =QFS.l(f);
     y2(h,draw,"cobb") = Y.l(h);

     display fd2, qp2, qc2, qx2, pm2, wt2, wk2, qfs2 ;


* 2) solve again with CES
*----------------------------------
* reset calibrated values
* --------------------
     FD.l(g,f,h)    = fd1(g,f,h,draw,"ces") ;
     QP.l(g,h)      = qp1(g,h,draw,"ces");
     QC.l(g,h)      = qc1(g,h,draw,"ces");
     QX.l(g)        = qx1(g,draw,"ces");
     WK.l(g,fk,h)   = wk1(g,fk,h,draw,"ces") ;
     WT.l(ft)       = wt1(ft,draw,"ces") ;
     PM.l(g)        = pm1(g,draw,"ces") ;
     QFS.l(ft)      = qfs1(ft,draw,"ces") ;

*and solve
     solve genCES using mcp ;
     abort$(genCES.modelstat ne 1) "ces simulation fails" ;
     display "post CES"  ;
     display FD.l, QP.l, QC.l, QX.l, PM.l, WT.l, WK.l, Y.l ;

     fd2(g,f,h,draw,"ces")    = FD.l(g,f,h);
     qp2(g,h,draw,"ces")      = QP.l(g,h);
     qc2(g,h,draw,"ces")      = QC.l(g,h);
     qx2(g,draw,"ces")        = QX.l(g);
     pm2(g,draw,"ces")        = PM.l(g);
     wt2(f,draw,"ces")        = WT.l(f);
     wk2(g,f,h,draw,"ces")    = WK.l(g,f,h);
     qfs2(f,draw,"ces")       = QFS.l(f);
     y2(h,draw,"ces")         = Y.l(h);
);

* output deltas and percentage changes parameters
fdD(g,f,h,draw,mod) = fd2(g,f,h,draw,mod)-fd1(g,f,h,draw,mod);
qpD(g,h,draw,mod)   = qp2(g,h,draw,mod)-qp1(g,h,draw,mod);
qcD(g,h,draw,mod)   = qc2(g,h,draw,mod)-qc1(g,h,draw,mod);
qxD(g,draw,mod)     = qx2(g,draw,mod)-qx1(g,draw,mod);
pmD(g,draw,mod)     = pm2(g,draw,mod)-pm1(g,draw,mod);
wtD(f,draw,mod)     = wt2(f,draw,mod)-wt1(f,draw,mod);
wkD(g,f,h,draw,mod) = wk2(g,f,h,draw,mod)-wk1(g,f,h,draw,mod);
qfsD(f,draw,mod)    = qfs2(f,draw,mod)-qfs1(f,draw,mod);
yD(h,draw,mod)      = y2(h,draw,mod)-y1(h,draw,mod);

fdPC(g,f,h,draw,mod)$fd1(g,f,h,draw,mod)     = 100* fdD(g,f,h,draw,mod)/fd1(g,f,h,draw,mod);
qpPC(g,h,draw,mod)$qp1(g,h,draw,mod)         = 100* qpD(g,h,draw,mod)/qp1(g,h,draw,mod);
qcPC(g,h,draw,mod)$qc1(g,h,draw,mod)         = 100* qcD(g,h,draw,mod)/qc1(g,h,draw,mod);
qxPC(g,draw,mod)$qx1(g,draw,mod)             = 100* qxD(g,draw,mod)/qx1(g,draw,mod);
pmPC(g,draw,mod)$pm1(g,draw,mod)             = 100* pmD(g,draw,mod)/pm1(g,draw,mod);
wtPC(f,draw,mod)$wt1(f,draw,mod)             = 100* wtD(f,draw,mod)/wt1(f,draw,mod);
wkPC(g,f,h,draw,mod)$wk1(g,f,h,draw,mod)     = 100* wkD(g,f,h,draw,mod)/wk1(g,f,h,draw,mod);
qfsPC(f,draw,mod)$qfs1(f,draw,mod)           = 100* qfsD(f,draw,mod)/qfs1(f,draw,mod);
yPC(h,draw,mod)$y1(h,draw,mod)               = 100* yD(h,draw,mod)/y1(h,draw,mod);

option fdPC:3:2:2;
option qpPC:3:1:2;
option qcPC:3:1:2;

display fdD, qpD, qcD, qxD, pmD, wtD, wkD, qfsD ;
display fdPC, qpPC, qcPC, qxPC, pmPC, wtPC, wkPC, qfsPC ;


*-----------------------------------------------------------------------
* Now output parameters with mean and variance
*-----------------------------------------------------------------------
set mv /mean, stdev/ ;

parameter
fdPC_mv(g,f,h,mv,mod)    mean and sample stdev of % change in fd over the draws.
qfsPC_mv(f,mv,mod)       mean and sample stdev of % change in qfs over the draws.
qpPC_mv(g,h,mv,mod)      mean and sample stdev of % change in qp over the draws.
qcPC_mv(g,h,mv,mod)      mean and sample stdev of % change in qc over the draws.
qxPC_mv(g,mv,mod)        mean and sample stdev of % change in qx over the draws.
pmPC_mv(g,mv,mod)        mean and sample stdev of % change in pm over the draws.
wtPC_mv(f,mv,mod)        mean and sample stdev of % change in wt over the draws.
wkPC_mv(g,f,h,mv,mod)    mean and sample stdev of % change in wk over the draws.
yPC_mv(h,mv,mod)         mean and sample stdev of % change in y over the draws.
;

fdPC_mv(g,f,h,"mean",mod) = sum(draw, fdPC(g,f,h,draw,mod)) / card(draw) ;
fdPC_mv(g,f,h,"stdev",mod) = sqrt(sum(draw, sqr(fdPC(g,f,h,draw,mod) - fdPC_mv(g,f,h,"mean",mod)))/(card(draw)-1)) ;

qfsPC_mv(f,"mean",mod) = sum(draw, qfsPC(f,draw,mod)) / card(draw) ;
qfsPC_mv(f,"stdev",mod) = sqrt(sum(draw, sqr(qfsPC(f,draw,mod) - qfsPC_mv(f,"mean",mod)))/(card(draw)-1)) ;

qpPC_mv(g,h,"mean",mod) = sum(draw, qpPC(g,h,draw,mod)) / card(draw) ;
qpPC_mv(g,h,"stdev",mod) = sqrt(sum(draw, sqr(qpPC(g,h,draw,mod) - qpPC_mv(g,h,"mean",mod)))/(card(draw)-1)) ;

qcPC_mv(g,h,"mean",mod) = sum(draw, qcPC(g,h,draw,mod)) / card(draw) ;
qcPC_mv(g,h,"stdev",mod) = sqrt(sum(draw, sqr(qcPC(g,h,draw,mod) - qcPC_mv(g,h,"mean",mod)))/(card(draw)-1)) ;

qxPC_mv(g,"mean",mod) = sum(draw, qxPC(g,draw,mod)) / card(draw) ;
qxPC_mv(g,"stdev",mod) = sqrt(sum(draw, sqr(qxPC(g,draw,mod) - qxPC_mv(g,"mean",mod)))/(card(draw)-1)) ;

pmPC_mv(g,"mean",mod) = sum(draw, pmPC(g,draw,mod)) / card(draw) ;
pmPC_mv(g,"stdev",mod) = sqrt(sum(draw, sqr(pmPC(g,draw,mod) - pmPC_mv(g,"mean",mod)))/(card(draw)-1)) ;

wkPC_mv(g,f,h,"mean",mod) = sum(draw, wkPC(g,f,h,draw,mod)) / card(draw) ;
wkPC_mv(g,f,h,"stdev",mod) = sqrt(sum(draw, sqr(wkPC(g,f,h,draw,mod) - wkPC_mv(g,f,h,"mean",mod)))/(card(draw)-1)) ;

wtPC_mv(f,"mean",mod) = sum(draw, wtPC(f,draw,mod)) / card(draw) ;
wtPC_mv(f,"stdev",mod) = sqrt(sum(draw, sqr(wtPC(f,draw,mod) - wtPC_mv(f,"mean",mod)))/(card(draw)-1)) ;

yPC_mv(h,"mean",mod) = sum(draw, yPC(h,draw,mod)) / card(draw) ;
yPC_mv(h,"stdev",mod) = sqrt(sum(draw, sqr(yPC(h,draw,mod) - yPC_mv(h,"mean",mod)))/(card(draw)-1)) ;


display fdPC_mv, qfsPC_mv, qpPC_mv, qcPC_mv, qxPC_mv, pmPC_mv, wtPC_mv, wkPC_mv, yPC_mv ;


