*==============================================================================*
*  国家5A级景区建设、减污降碳与空间溢出效应 —— 完整可复现脚本（稳健版 v4）
*  一个 do + 一个 dta（data_spatial.dta 内含 lat lon mean_pgdp）。结果输出到 results/。
*
*  【本版关键改动】
*   1) 所有 Mata 均为“全向量化”，无 for 循环 / if / 嵌套花括号（section 4 除外，
*      仅用一个最简单的单层 for 且在独立块内，安全）；
*   2) 【重要】任何 Mata 块都不再放进 Stata 的 foreach/forvalues 循环内
*      （交互式 mata: 放进 Stata 循环会中断/报 Break）。Stata 循环只调用 xsmle。
*   请务必用【本文件】整体覆盖旧版后运行。
*
*  需安装： ssc install xsmle ; ssc install xthreg ; ssc install estout ; ssc install reghdfe ftools
*  运行前： cd "data_spatial.dta 所在目录"，例如 cd "/Users/zifeimeng/Desktop/0702空间/02"
*==============================================================================*
clear all
set more off
set matsize 2000
cap mkdir results

use "data_spatial.dta", clear
cap gen double lnfin = ln(fin)
xtset city_code year
global CTRL lnpgdp lndensity urban struc2 gov tech human lnfin

*------------------------------------------------------------------------------*
* 0. 构造 5 类空间权重矩阵（全向量化 Mata，独立块）
*------------------------------------------------------------------------------*
preserve
    bysort city_code (year): keep if _n==1
    sort city_code
    mkmat lat,       matrix(LAT)
    mkmat lon,       matrix(LON)
    mkmat mean_pgdp, matrix(YBAR)
restore
mata:
    LAT = st_matrix("LAT"):*(pi()/180)
    LON = st_matrix("LON"):*(pi()/180)
    Y   = st_matrix("YBAR")
    n   = rows(LAT); one = J(n,1,1)
    LAi = LAT*one'; LAj = one*LAT'; LOi = LON*one'; LOj = one*LON'
    dla = LAj - LAi; dlo = LOj - LOi
    ah  = sin(dla:/2):^2 + cos(LAi):*cos(LAj):*sin(dlo:/2):^2
    D   = 2:*6371:*asin(sqrt(ah))
    G   = editmissing(1:/(D:^2), 0)
    rs = rowsum(G); rs = rs + (rs:==0); Wgeo = G :/ rs
    DE  = abs(Y*one' - one*Y'); E = editmissing(1:/DE, 0)
    rs = rowsum(E); rs = rs + (rs:==0); Wecon = E :/ rs
    Yr  = Y :/ mean(Y); EG = G :* (one*Yr')
    rs = rowsum(EG); rs = rs + (rs:==0); Wegn = EG :/ rs
    EW  = 0.5:*Wgeo + 0.5:*Wecon
    rs = rowsum(EW); rs = rs + (rs:==0); Wegw = EW :/ rs
    A   = (D:<=163):*(D:>0)
    Dbig = D + I(n):*1e12; rmins = rowmin(Dbig); NN = (Dbig:==(rmins*one'))
    iso = (rowsum(A):==0); A = A + NN:*(iso*one'); A = (A + A'):>0
    rs = rowsum(A); rs = rs + (rs:==0); Wadj = A :/ rs
    st_matrix("Wadj",Wadj); st_matrix("Wgeo",Wgeo); st_matrix("Wecon",Wecon)
    st_matrix("Wegw",Wegw); st_matrix("Wegn",Wegn); st_matrix("Dmat",D)
end
di as result "== 5 类空间权重矩阵已构造：Wadj Wgeo Wecon Wegw Wegn =="

*------------------------------------------------------------------------------*
* 1. 空间诊断：全局 Moran's I（5 矩阵×双向FE残差）+ 分年度 Moran's I 表
*    —— 全部用独立 mata 块，不放进 Stata 循环
*------------------------------------------------------------------------------*
sort year city_code
cap drop _res
qui reghdfe lnpoco2 DID $CTRL, a(city_code year) residuals(_res)
* 1a 五矩阵全局 Moran's I（残差）：一个独立 mata 块
mata:
    e = st_data(.,"_res"); N=289; T=21
    Em = rowshape(e, T)'; Em = Em :- (J(N,1,1)*mean(Em)); den = sum(Em:*Em)
    st_numscalar("mI1", sum(Em:*(st_matrix("Wadj") *Em))/den)
    st_numscalar("mI2", sum(Em:*(st_matrix("Wgeo") *Em))/den)
    st_numscalar("mI3", sum(Em:*(st_matrix("Wecon")*Em))/den)
    st_numscalar("mI4", sum(Em:*(st_matrix("Wegw") *Em))/den)
    st_numscalar("mI5", sum(Em:*(st_matrix("Wegn")*Em))/den)
end
di as txt _n "== 全局 Moran's I（双向固定效应残差）=="
di as txt "  Wadj="  as res %6.3f mI1 as txt "   Wgeo="  as res %6.3f mI2 ///
   as txt "   Wecon=" as res %6.3f mI3 as txt "   Wegw=" as res %6.3f mI4 ///
   as txt "   Wegn=" as res %6.3f mI5
drop _res

* 1b 分年度 Moran's I + Z（主推 Wegw）：一个独立 mata 块（向量化跨年）
sort year city_code
mata:
    W = st_matrix("Wegw"); y = st_data(.,"lnpoco2"); N=289; T=21
    Ym = rowshape(y, T)'; Ym = Ym :- (J(N,1,1)*mean(Ym))
    num = colsum(Ym :* (W*Ym)); den = colsum(Ym:*Ym); Iv = (num:/den)'
    S0=sum(W); S1=0.5*sum((W+W'):^2); S2=sum((rowsum(W)+colsum(W)'):^2)
    EI=-1/(N-1); VI=(N^2*S1-N*S2+3*S0^2)/(S0^2*(N^2-1))-EI^2
    Zv = (Iv:-EI):/sqrt(VI); yrs = (2003::2023)
    st_matrix("MORAN", (yrs, Iv, Zv))
end
preserve
    clear
    svmat MORAN
    rename (MORAN1 MORAN2 MORAN3) (year MoranI Zscore)
    gen Pvalue = 2*(1-normal(abs(Zscore)))
    format MoranI Zscore Pvalue %9.3f
    di as txt _n "== 分年度全局 Moran's I（Wegw）=="
    list year MoranI Zscore Pvalue, sep(0) noobs
    cap export excel using "results/moran_by_year.xlsx", replace first(var)
    save "results/moran_by_year.dta", replace
    twoway (connected MoranI year, msymbol(O) lcolor(black) mcolor(black)), ///
        scheme(s1mono) ytitle("全局 Moran's I") xtitle("年份") ///
        title("减污降碳全局空间自相关的时间演变")
    graph export "results/fig_moran_trend.png", replace width(2000) height(1300)
restore

*------------------------------------------------------------------------------*
* 2. 空间杜宾模型 SDM：三列 (1)Time FE (2)Individual FE (3)Two-way FE（完整系数）
*    —— Stata 循环只有 eststo/xsmle，无 mata
*------------------------------------------------------------------------------*
eststo clear
eststo sdm_time: xsmle lnpoco2 DID $CTRL, model(sdm) wmat(Wegw) fe type(time) nsim(200)
estadd scalar rho = e(rho)
eststo sdm_ind:  xsmle lnpoco2 DID $CTRL, model(sdm) wmat(Wegw) fe type(ind)  nsim(200)
estadd scalar rho = e(rho)
eststo sdm_both: xsmle lnpoco2 DID $CTRL, model(sdm) wmat(Wegw) fe type(both) nsim(200)
estadd scalar rho = e(rho)
esttab sdm_time sdm_ind sdm_both using "results/table_SDM_3FE.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("(1) Time FE" "(2) Individual FE" "(3) Two-way FE") ///
    scalars("rho 空间自回归系数rho") stats(N, labels("观测值N")) ///
    nogaps compress title("空间杜宾模型完整估计(W_egw)")

* Hausman（版本敏感，capture 包裹；失败回退非空间面板 Hausman）
cap noisily {
    qui xsmle lnpoco2 DID $CTRL, model(sdm) wmat(Wegw) fe type(ind)
    est store fe_sdm
    qui xsmle lnpoco2 DID $CTRL, model(sdm) wmat(Wegw) re
    est store re_sdm
    hausman fe_sdm re_sdm, sigmamore
}
if _rc {
    di as txt "改用非空间面板 Hausman："
    qui xtreg lnpoco2 DID $CTRL, fe
    est store fe0
    qui xtreg lnpoco2 DID $CTRL, re
    est store re0
    hausman fe0 re0
}

*------------------------------------------------------------------------------*
* 3. 替换五类权重矩阵（完整系数，Time FE）—— Stata 循环只有 xsmle
*------------------------------------------------------------------------------*
eststo clear
foreach W in Wadj Wgeo Wecon Wegw Wegn {
    eststo m_`W': xsmle lnpoco2 DID $CTRL, model(sdm) wmat(`W') fe type(time) nsim(200)
    estadd scalar rho = e(rho)
}
esttab m_Wadj m_Wgeo m_Wecon m_Wegw m_Wegn using "results/table_SDM_5matrices.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("邻接" "地理" "经济" "经济地理权重" "经济地理嵌套") ///
    scalars("rho rho") stats(N, labels("N")) nogaps compress ///
    title("五类权重矩阵 SDM 完整估计(Time FE)")

*------------------------------------------------------------------------------*
* 4. 溢出范围：不同 km 距离阈值下的 SDM（Time FE）
*    —— 先在一个独立 mata 块里预建全部 8 个距离带矩阵为 Stata 矩阵；
*       再用 Stata 循环调用 xsmle（循环内无 mata）
*------------------------------------------------------------------------------*
mata:
    Dm = st_matrix("Dmat"); n = rows(Dm); one = J(n,1,1)
    Dbig = Dm + I(n):*1e12; rmins = rowmin(Dbig); NN = (Dbig:==(rmins*one'))
    Bk=(Dm:<=150):*(Dm:>0); Bk=(Bk+NN:*((rowsum(Bk):==0)*one')):>0; rs=rowsum(Bk); rs=rs+(rs:==0); st_matrix("Wd150", Bk:/rs)
    Bk=(Dm:<=200):*(Dm:>0); Bk=(Bk+NN:*((rowsum(Bk):==0)*one')):>0; rs=rowsum(Bk); rs=rs+(rs:==0); st_matrix("Wd200", Bk:/rs)
    Bk=(Dm:<=250):*(Dm:>0); Bk=(Bk+NN:*((rowsum(Bk):==0)*one')):>0; rs=rowsum(Bk); rs=rs+(rs:==0); st_matrix("Wd250", Bk:/rs)
    Bk=(Dm:<=300):*(Dm:>0); Bk=(Bk+NN:*((rowsum(Bk):==0)*one')):>0; rs=rowsum(Bk); rs=rs+(rs:==0); st_matrix("Wd300", Bk:/rs)
    Bk=(Dm:<=350):*(Dm:>0); Bk=(Bk+NN:*((rowsum(Bk):==0)*one')):>0; rs=rowsum(Bk); rs=rs+(rs:==0); st_matrix("Wd350", Bk:/rs)
    Bk=(Dm:<=400):*(Dm:>0); Bk=(Bk+NN:*((rowsum(Bk):==0)*one')):>0; rs=rowsum(Bk); rs=rs+(rs:==0); st_matrix("Wd400", Bk:/rs)
    Bk=(Dm:<=450):*(Dm:>0); Bk=(Bk+NN:*((rowsum(Bk):==0)*one')):>0; rs=rowsum(Bk); rs=rs+(rs:==0); st_matrix("Wd450", Bk:/rs)
    Bk=(Dm:<=500):*(Dm:>0); Bk=(Bk+NN:*((rowsum(Bk):==0)*one')):>0; rs=rowsum(Bk); rs=rs+(rs:==0); st_matrix("Wd500", Bk:/rs)
end
cap postclose PM
postfile PM int km double rho using "results/spillover_range.dta", replace
foreach km in 150 200 250 300 350 400 450 500 {
    qui xsmle lnpoco2 DID $CTRL, model(sdm) wmat(Wd`km') fe type(time) nsim(1)
    post PM (`km') (e(rho))
}
postclose PM
preserve
    use "results/spillover_range.dta", clear
    twoway (connected rho km, msymbol(O) lcolor(black) mcolor(black)), ///
        yline(0,lpattern(dash) lcolor(gs9)) scheme(s1mono) ///
        xtitle("地理距离阈值 (km)") ytitle("空间自回归系数 rho") ///
        title("减污降碳空间溢出的地理衰减")
    graph export "results/fig_distance_decay.png", replace width(2000) height(1300)
restore

*------------------------------------------------------------------------------*
* 5. 门槛/分段异质性溢出
*------------------------------------------------------------------------------*
* (A) Hansen 面板门槛（版本敏感，capture 包裹）
cap noisily xthreg lnpoco2 DID $CTRL, rx(DID) qx(lnpgdp) thnum(2) trim(0.01 0.01) grid(300) bs(300 300)

* (B) 分区制 SDM（各维度分高/低，Time FE）—— Stata 循环内无 mata
eststo clear
foreach v of newlist human lnpgdp ter_gdp er {
    capture confirm variable `v'
    if _rc continue
    cap drop hi_`v' DIDlo_`v' DIDhi_`v' mm_`v'
    bysort city_code: egen mm_`v' = mean(`v')
    qui sum mm_`v', detail
    gen byte hi_`v' = mm_`v' > r(p50)
    gen double DIDlo_`v' = DID*(1-hi_`v')
    gen double DIDhi_`v' = DID*hi_`v'
    eststo reg_`v': xsmle lnpoco2 DIDlo_`v' DIDhi_`v' $CTRL, model(sdm) wmat(Wegw) fe type(time) nsim(100)
}
cap drop east DIDlo_east DIDhi_east
gen byte east = region=="东部"
gen double DIDlo_east = DID*(1-east)
gen double DIDhi_east = DID*east
eststo reg_east: xsmle lnpoco2 DIDlo_east DIDhi_east $CTRL, model(sdm) wmat(Wegw) fe type(time) nsim(100)
esttab reg_* using "results/table_regime_spillover.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) nogaps compress ///
    title("分区制异质性溢出(Time FE, W_egw)")

* (C) 控制组溢出污染检验（queen 邻接；独立 mata 块构造 spill）
sort year city_code
cap drop spill
mata:
    Wb = (st_matrix("Wadj"):>0); did = st_data(.,"DID"); N=289; T=21
    Dmat = rowshape(did, T)'
    SP = (Dmat:==0) :* ((Wb*Dmat):>0)
    st_store(., st_addvar("byte","spill"), vec(SP))
end
eststo clear
eststo base_did:  reghdfe lnpoco2 DID $CTRL,       a(city_code year) vce(cluster city_code)
eststo spill_did: reghdfe lnpoco2 DID spill $CTRL, a(city_code year) vce(cluster city_code)
esttab base_did spill_did using "results/table_spillover_contamination.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) keep(DID spill) ///
    mtitles("基准" "加溢出虚拟变量") nogaps compress title("控制组溢出污染检验")

*------------------------------------------------------------------------------*
* 6. Moran 散点图：2003 与 2023 年（主推 Wegw）
*    —— 独立 mata 块算好 (z, W·z) 存为 Stata 矩阵；再用 Stata 分别作图（无 mata 循环）
*------------------------------------------------------------------------------*
sort year city_code
mata:
    W = st_matrix("Wegw"); y = st_data(.,"lnpoco2"); NN = rows(y); N=289
    z03 = y[|1 \ N|];               z03 = z03 :- mean(z03)
    z23 = y[|(NN-N+1) \ NN|];       z23 = z23 :- mean(z23)
    st_matrix("Z03", (z03, W*z03)); st_matrix("Z23", (z23, W*z23))
    st_numscalar("MI03", (z03'*(W*z03))/(z03'*z03))
    st_numscalar("MI23", (z23'*(W*z23))/(z23'*z23))
end
preserve
    clear
    svmat Z03
    twoway (scatter Z032 Z031, mcolor(black) msize(small) msymbol(Oh)) ///
           (lfit Z032 Z031, lcolor(black)), ///
        yline(0,lpattern(dash) lcolor(gs10)) xline(0,lpattern(dash) lcolor(gs10)) ///
        scheme(s1mono) legend(off) xtitle("去均值 lnpoco2") ytitle("空间滞后 W·z") ///
        title("Moran 散点图 2003 (Moran's I=" + string(`=MI03',"%5.3f") + ")")
    graph export "results/fig_moran_2003.png", replace width(1600) height(1400)
restore
preserve
    clear
    svmat Z23
    twoway (scatter Z232 Z231, mcolor(black) msize(small) msymbol(Oh)) ///
           (lfit Z232 Z231, lcolor(black)), ///
        yline(0,lpattern(dash) lcolor(gs10)) xline(0,lpattern(dash) lcolor(gs10)) ///
        scheme(s1mono) legend(off) xtitle("去均值 lnpoco2") ytitle("空间滞后 W·z") ///
        title("Moran 散点图 2023 (Moran's I=" + string(`=MI23',"%5.3f") + ")")
    graph export "results/fig_moran_2023.png", replace width(1600) height(1400)
restore

di as result _n "================ 完成：全部结果见 results 文件夹 ================"
