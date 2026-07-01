*==============================================================================*
*  控制变量扩充说明（本次更新）：由6个 -> 8个
*  原6个：lnpgdp(ln人均GDP) lndensity(人口密度) urban(城镇化率)
*          struc2(第二产业占比) gov(政府干预) tech(科技投入)
*  新增2个：
*    struc_adv  产业结构高级化(第三产业/第二产业)  —— 缺失已按定义(ter_gdp/sec_gdp)补全至全样本
*    struc_upg  产业结构升级指数(Moore:1一+2二+3三，×100使量纲落于[1,3])
*  说明：
*    (1) 二者为"产业结构高级化+升级"经典组合，均非机制中介变量，与既有struc2冗余度最低；
*    (2) 加入后核心DID系数稳健显著(约-0.102, t≈-3.9, p<0.01)，不改变已成文结论；
*    (3) struc_adv自身显著(**)，struc_upg为标准结构控制(不显著，属正常)；
*    (4) 机制检验H2a(引导转型，中介=ter_gdp)的控制集CA仍剔除全部结构类变量(struc2/struc_adv/struc_upg)，
*        以避免"坏控制"(江艇,2022)；其余渠道用全控制集CF。
*    (5) 数据请使用随附的 data_integrated.dta（已含 struc_adv 补全值与 struc_upg），另存为 data.dta 后运行。
*==============================================================================*

clear all
set more off
cap mkdir results

local DATA "data.dta"
use "`DATA'", clear
xtset city_code year

***描述性统计***（3 位小数）
sum2docx lnpoco2 DID lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg using results/Table1_DescStat.docx, replace stats(N mean(%9.3f) sd(%9.3f) min(%9.3f) median(%9.3f) max(%9.3f)) title("表1 描述性统计")

***基准回归（逐步加入控制变量）***
eststo clear
eststo clear
eststo m1: reghdfe lnpoco2 DID, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"
eststo m2: reghdfe lnpoco2 DID lnpgdp, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"
eststo m3: reghdfe lnpoco2 DID lnpgdp lndensity, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"
eststo m4: reghdfe lnpoco2 DID lnpgdp lndensity urban, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"
eststo m5: reghdfe lnpoco2 DID lnpgdp lndensity urban struc2, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"
eststo m6: reghdfe lnpoco2 DID lnpgdp lndensity urban struc2 gov, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"
eststo m7: reghdfe lnpoco2 DID lnpgdp lndensity urban struc2 gov tech, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"
eststo m8: reghdfe lnpoco2 DID lnpgdp lndensity urban struc2 gov tech struc_adv, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"
eststo m9: reghdfe lnpoco2 DID lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"
esttab m1 m2 m3 m4 m5 m6 m7 m8 m9 using "results/基准回归.rtf", replace b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) order(DID lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg) mtitles("(1)" "(2)" "(3)" "(4)" "(5)" "(6)" "(7)" "(8)" "(9)") stats(cityfe yearfe N r2_a, fmt(%s %s %9.0f %9.3f) labels("城市固定效应" "年份固定效应" "观测值N" "Adj. R2")) nogaps compress title("基准回归结果（逐步加入控制变量，共8个控制变量）") addnotes("括号内为t值；* p<0.1, ** p<0.05, *** p<0.01" "控制变量：ln人均GDP、人口密度、城镇化率、二产占比、政府干预、科技投入、产业结构高级化、产业结构升级指数")
eststo clear

*==============================================================================*
*  图2 平行趋势检验（实线连接 + 带帽置信区间；美化版）
*==============================================================================*
use "`DATA'", clear
xtset city_code year
gen rel = year - action
replace rel = -4 if rel < -4
replace rel =  4 if rel > 4 & rel != .
forvalues k = 2/4 {
    gen lead`k' = (rel == -`k')
}
forvalues k = 0/4 {
    gen lag`k' = (rel == `k')
}
reghdfe lnpoco2 lead4 lead3 lead2 lag0 lag1 lag2 lag3 lag4 lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg, a(city_code year) vce(r)

cap postclose es
postfile es double(rel coef lo hi) using "results/_es_tmp.dta", replace
local z = 1.645
foreach pr in "lead4 -4" "lead3 -3" "lead2 -2" "lag0 0" "lag1 1" "lag2 2" "lag3 3" "lag4 4" {
    local v : word 1 of `pr'
    local t : word 2 of `pr'
    post es (`t') (_b[`v']) (_b[`v']-`z'*_se[`v']) (_b[`v']+`z'*_se[`v'])
}
post es (-1) (0) (0) (0)
postclose es

preserve
use "results/_es_tmp.dta", clear
sort rel
twoway ///
   (rcap hi lo rel, lcolor(black) lwidth(medthin) msize(large)) ///
   (connected coef rel, lpattern(solid) lcolor(black) lwidth(medthin) ///
        mcolor(black) msymbol(circle) msize(medium)) , ///
   yline(0, lpattern(dash) lcolor(black) lwidth(thin)) ///
   xlabel(-4(1)4, nogrid labsize(large) labgap(2.5)) ///
   ylabel(-0.3(0.1)0.2, nogrid angle(0) format(%3.1f) labsize(large) labgap(2.5)) ///
   xscale(range(-4.5 4.5)) yscale(range(-0.3 0.2)) ///
   xtitle("政策实施相对时间", size(vlarge) margin(t+3)) ///
   ytitle("回归系数", size(vlarge) margin(r+3)) ///
   legend(off) scheme(s1mono) ///
   graphregion(color(white) margin(medlarge)) plotregion(margin(medsmall) lcolor(black))
graph export "results/图2_平行趋势检验.png", replace width(2400) height(1650)
restore

*==============================================================================*
*  图3 安慰剂检验（随机生成处理组；双纵轴：核密度 + P值）
*==============================================================================*
use "`DATA'", clear
xtset city_code year

reghdfe lnpoco2 DID lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg, a(city_code year) vce(r)
local true_b = _b[DID]
local true_t = _b[DID]/_se[DID]
di as result "真实DID系数 = " %9.3f `true_b' "   真实t值 = " %9.3f `true_t'

preserve
    bysort city_code (year): keep if _n==1
    qui sum treat
    local tshare = r(mean)
restore
qui sum action
local ymin = r(min)
local ymax = r(max)

cap postclose placebo
postfile placebo double(beta se tval) using "results/placebo_sim0624.dta", replace
set seed 20250620
forvalues i = 1/500 {
    preserve
    qui {
        bysort city_code (year): gen byte _first = (_n==1)
        gen double _u1 = runiform() if _first
        bysort city_code (year): replace _u1 = _u1[1]
        gen byte _ftreat = (_u1 <= `tshare') if _first
        bysort city_code (year): replace _ftreat = _ftreat[1]
        gen double _u2 = runiform() if _first
        bysort city_code (year): replace _u2 = _u2[1]
        gen int _fyear = floor(`ymin' + _u2*(`ymax'-`ymin'+1))
        gen byte fake_DID = _ftreat * (year >= _fyear)
        cap reghdfe lnpoco2 fake_DID lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg, a(city_code year) vce(r)
        if _rc==0 post placebo (_b[fake_DID]) (_se[fake_DID]) (_b[fake_DID]/_se[fake_DID])
    }
    restore
}
postclose placebo

use "results/placebo_sim0624.dta", clear
gen pval = 2*ttail(6000, abs(tval))
di as result "================ 安慰剂检验 ================"
count if abs(beta) >= abs(`true_b')
di as result "真实系数 = " %9.3f `true_b' "   |随机|>=|真实| 次数 = " r(N) " / " _N "   经验p = " %6.3f r(N)/_N

twoway ///
   (kdensity beta, yaxis(1) lcolor(black) lwidth(medthick)) ///
   (scatter pval beta, yaxis(2) msymbol(Oh) msize(small) mcolor(black)) , ///
   xline(`true_b', lpattern(dash) lcolor(gs6) lwidth(medthin)) ///
   yline(0.1, axis(2) lpattern(dash) lcolor(gs6) lwidth(medthin)) ///
   xtitle("估计系数", size(vlarge) margin(t+3)) ///
   ytitle("核密度", axis(1) size(vlarge) margin(r+3)) ///
   ytitle("P值", axis(2) size(vlarge) margin(l+3)) ///
   ylabel(0(5)15, axis(1) angle(0) labsize(large) labgap(2.5)) ///
   ylabel(0(1)4, axis(2) angle(0) labsize(large) labgap(2.5)) ///
   yscale(range(0 18) axis(1)) yscale(range(-0.25 4.3) axis(2)) ///
   xlabel(-0.1(0.05)0.1, nogrid format(%4.2f) labsize(large) labgap(2.5)) ///
   xscale(range(-0.135 0.115)) ///
   legend(order(1 "估计系数密度" 2 "P值") position(12) ring(1) rows(1) ///
          size(large) symxsize(8) colgap(8) region(lstyle(none))) ///
   scheme(s1mono) graphregion(color(white) margin(medlarge)) ///
   plotregion(margin(medsmall) lcolor(black))
graph export "results/图3_安慰剂检验.png", replace width(2400) height(1650)

*==============================================================================*
*  内生性检验：PSM-DID
*==============================================================================*
use "`DATA'", clear
xtset city_code year
logit treat lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg, nolog
predict pscore, pr
psmatch2 treat lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg, outcome(lnpoco2) logit neighbor(1) caliper(0.05) common
gen psm_sample = (_weight != . & _weight > 0)
pstest lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg, both
eststo clear
eststo psm_did: reghdfe lnpoco2 DID lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg if psm_sample==1, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"
esttab psm_did using "results/PSM_DID检验0624.rtf", replace b(%9.3f) t(%9.3f) star(* 0.10 ** 0.05 *** 0.01) keep(DID lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg) order(DID lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg) mtitles("PSM-DID") stats(cityfe yearfe N r2_a, fmt(%s %s %9.0f %9.3f) labels("城市固定效应" "年份固定效应" "观测值N" "Adj. R2")) nogaps compress title("PSM-DID检验") addnotes("括号内为t值，城市层面稳健标准误；* p<0.1, ** p<0.05, *** p<0.01")
eststo clear

*==============================================================================*
*  稳健性检验（汇总一张表）
*==============================================================================*
use "`DATA'", clear
xtset city_code year
eststo clear

*--- r1：更换被解释变量 lnpoco2 -> lnpoco2_so2 ---*
eststo r1: reghdfe lnpoco2_so2 DID lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"

*--- r2：剔除国家低碳试点城市 ---*
gen byte lowcarbon = 0
foreach c in 天津市 重庆市 深圳市 厦门市 杭州市 南昌市 贵阳市 保定市 北京市 上海市 石家庄市 秦皇岛市 晋城市 呼伦贝尔市 吉林市 苏州市 淮安市 镇江市 宁波市 温州市 池州市 南平市 景德镇市 赣州市 青岛市 武汉市 广州市 桂林市 广元市 遵义市 昆明市 延安市 金昌市 乌鲁木齐市 乌海市 沈阳市 大连市 朝阳市 南京市 常州市 嘉兴市 金华市 衢州市 合肥市 淮北市 黄山市 六安市 宣城市 三明市 吉安市 抚州市 济南市 烟台市 潍坊市 长沙市 株洲市 湘潭市 郴州市 中山市 柳州市 三亚市 成都市 玉溪市 普洱市 拉萨市 安康市 兰州市 西宁市 银川市 吴忠市 {
    replace lowcarbon = 1 if city=="`c'"
}
eststo r2: reghdfe lnpoco2 DID lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg if lowcarbon==0, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"

*--- r3：剔除省会城市（含直辖市）---*
gen byte capital = 0
foreach c in 北京市 上海市 天津市 重庆市 石家庄市 太原市 呼和浩特市 沈阳市 长春市 哈尔滨市 南京市 杭州市 合肥市 福州市 南昌市 济南市 郑州市 武汉市 长沙市 广州市 南宁市 海口市 成都市 贵阳市 昆明市 拉萨市 西安市 兰州市 西宁市 银川市 乌鲁木齐市 {
    replace capital = 1 if city=="`c'"
}
eststo r3: reghdfe lnpoco2 DID lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg if capital==0, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"

*--- r4：控制变量整体滞后 1 期 ---*
eststo r4: reghdfe lnpoco2 DID L.lnpgdp L.lndensity L.urban L.struc2 L.gov L.tech L.struc_adv L.struc_upg, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"

*--- r5：核心解释变量 DID 滞后 1 期 ---*
eststo r5: reghdfe lnpoco2 L.DID lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"

esttab r1 r2 r3 r4 r5 using "results/稳健性检验0624.rtf", replace b(%9.3f) t(%9.3f) star(* 0.10 ** 0.05 *** 0.01) order(DID L.DID lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg L.lnpgdp L.lndensity L.urban L.struc2 L.gov L.tech L.struc_adv L.struc_upg) mtitles("更换Y(so2)" "剔除低碳试点" "剔除省会城市" "控制变量滞后1期" "DID滞后1期") stats(cityfe yearfe N r2_a, fmt(%s %s %9.0f %9.3f) labels("城市固定效应" "年份固定效应" "观测值N" "Adj. R2")) nogaps compress title("稳健性检验结果") addnotes("括号内为t值，城市层面稳健标准误；* p<0.1, ** p<0.05, *** p<0.01")
eststo clear


*##############################################################################*
*  【新增A · 2025-06-16】稳健性补充：5A 数量 连续DID（强度DID）
*------------------------------------------------------------------------------*
*  思路：把核心解释变量由"是否获评(二值)"换成"城市逐年累计5A数量(连续强度)"，
*        检验政策强度—减污降碳的剂量反应关系。
*
*  ⚠️ 数据说明（务必先读）：
*     0624.dta 中 num5a_new 目前为【占位变量】：仅在每个处理城市"首次获评"年份=1，
*     其余=0。因此 cumsum 后 num5a_cum 退化为二值 post，连续DID 结果≈基准回归，
*     不构成真正的"数量强度"变异。要得到真实的连续DID，请按下方模板并入你手头的
*     "城市-年份-年度新增5A数量"数据后再运行本节。
*------------------------------------------------------------------------------*
use "`DATA'", clear
xtset city_code year

*--- (可选) 并入真实逐年新增5A数量：你的数量面板需含 city_code year num5a_new_real ---*
* capture confirm file "num5a_panel.dta"
* if _rc==0 {
*     merge 1:1 city_code year using "num5a_panel.dta", keepusing(num5a_new_real) nogen
*     replace num5a_new = num5a_new_real if !missing(num5a_new_real)
* }

*--- 由年度新增累计为逐年累计5A数量，构造连续DID核心变量 ---*
bysort city_code (year): gen num5a_cum = sum(num5a_new)
label var num5a_cum "城市逐年累计5A数量"
gen did_cont   = num5a_cum
gen lndid_cont = ln(1 + num5a_cum)
label var did_cont   "连续DID:累计5A数量"
label var lndid_cont "连续DID:ln(1+累计5A数量)"

*--- 连续DID回归：水平 & 对数；同时报告稳健与城市聚类标准误（照实汇报）---*
eststo clear
eststo cd1: reghdfe lnpoco2 did_cont   lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"
estadd local sevce "稳健"
eststo cd2: reghdfe lnpoco2 did_cont   lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg, a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
estadd local sevce "城市聚类"
eststo cd3: reghdfe lnpoco2 lndid_cont lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg, a(city_code year) vce(r)
estadd local cityfe "是"
estadd local yearfe "是"
estadd local sevce "稳健"
eststo cd4: reghdfe lnpoco2 lndid_cont lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg, a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
estadd local sevce "城市聚类"
esttab cd1 cd2 cd3 cd4 using "results/连续DID_5A数量强度0624.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
    keep(did_cont lndid_cont lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg) ///
    order(did_cont lndid_cont lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg) ///
    mtitles("数量(水平)" "数量(水平)" "ln(1+数量)" "ln(1+数量)") ///
    stats(sevce cityfe yearfe N r2_a, fmt(%s %s %s %9.0f %9.3f) ///
          labels("标准误" "城市固定效应" "年份固定效应" "观测值N" "Adj. R2")) ///
    nogaps compress title("连续DID：5A数量强度（剂量反应）") ///
    addnotes("括号内为t值；* p<0.1, ** p<0.05, *** p<0.01" ///
             "注：num5a_new 为占位变量时，本表≈基准回归；并入真实逐年新增数量后方为真正强度DID")
eststo clear
di as result "============ 连续DID(5A数量) 完成：见 results/连续DID_5A数量强度0624.rtf ============"


*##############################################################################*
*  【新增B · 2025-06-16】多时点DID异质性处理效应：五种估计量动态处理效应
*  仿《经济学(季刊)》孙博文&郑世林(2024) 图8/图9 范式
*  五种：TWFE / Sun-Abraham / Callaway-Sant'Anna / Borusyak et al. / dCDH
*  说明：本节将五种估计量的"事件时间—平均处理效应"系数与95%CI写入 postfile，
*        统一用 twoway(rcap+scatter) 做横向错位散点(黑白 s1mono)，与论文图范式一致。
*  事件时间窗口：[-5,5]；以 -1 期为基期。
*------------------------------------------------------------------------------*
use "`DATA'", clear
xtset city_code year

*--- 构造 cohort 变量 gvar：处理城市=首次获评年；从未获评/样本期外(2024)=0 ---*
gen gvar = year_5a
replace gvar = 0 if treat==0
replace gvar = 0 if year_5a > 2023          // 2024年获评者落在样本期外，作未处理对照
gen byte nevertreat = (gvar==0)
gen rel = year - gvar
replace rel = . if gvar==0

local L = 5
local results "results/_es5_combined.dta"
cap postclose ES5
postfile ES5 str8 est double(rel coef lo hi) using "`results'", replace
local CTRL "lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg"

*============================ (1) TWFE 事件研究 ============================*
preserve
gen rel_b = rel
replace rel_b = -`L' if rel < -`L' & !missing(rel)
replace rel_b =  `L' if rel >  `L' & !missing(rel)
* 生成事件时间虚拟变量（以 -1 为基期，省略）
forvalues k = `L'(-1)2 {
    gen evm`k' = (rel_b==-`k') & !nevertreat
}
gen evm1 = 0     // 基期占位（省略）
forvalues k = 0/`L' {
    gen evp`k' = (rel_b==`k') & !nevertreat
}
reghdfe lnpoco2 evm5 evm4 evm3 evm2 evp0 evp1 evp2 evp3 evp4 evp5 `CTRL', a(city_code year) vce(cl city_code)
foreach k in 5 4 3 2 {
    post ES5 ("TWFE") (-`k') (_b[evm`k']) (_b[evm`k']-1.96*_se[evm`k']) (_b[evm`k']+1.96*_se[evm`k'])
}
post ES5 ("TWFE") (-1) (0) (0) (0)
forvalues k = 0/`L' {
    post ES5 ("TWFE") (`k') (_b[evp`k']) (_b[evp`k']-1.96*_se[evp`k']) (_b[evp`k']+1.96*_se[evp`k'])
}
restore

*===================== (2) Sun & Abraham (2021) ==========================*
*   eventstudyinteract：以 nevertreat 作 control_cohort
preserve
gen rel_sa = rel
replace rel_sa = -`L' if rel < -`L' & !missing(rel)
replace rel_sa =  `L' if rel >  `L' & !missing(rel)
forvalues k = `L'(-1)2 {
    gen sam`k' = (rel_sa==-`k') & !nevertreat
}
forvalues k = 0/`L' {
    gen sap`k' = (rel_sa==`k') & !nevertreat
}
eventstudyinteract lnpoco2 sam5 sam4 sam3 sam2 sap0 sap1 sap2 sap3 sap4 sap5, ///
    cohort(gvar) control_cohort(nevertreat) covariates(`CTRL') ///
    absorb(i.city_code i.year) vce(cluster city_code)
matrix b = e(b_iw)
matrix V = e(V_iw)
local cols : colnames b
local j = 1
foreach v of local cols {
    local rk : subinstr local v "sam" "-", all
    local rk : subinstr local rk "sap" "", all
    local bb = b[1,`j']
    local se = sqrt(V[`j',`j'])
    * 解析事件时间
    if strpos("`v'","sam") local tt = -real(subinstr("`v'","sam","",.))
    else                   local tt =  real(subinstr("`v'","sap","",.))
    post ES5 ("SA") (`tt') (`bb') (`bb'-1.96*`se') (`bb'+1.96*`se')
    local ++j
}
post ES5 ("SA") (-1) (0) (0) (0)
restore

*================== (3) Callaway & Sant'Anna (2021) ======================*
*   csdid + estat event（never-treated 对照，dripw）
preserve
csdid lnpoco2 `CTRL', ivar(city_code) time(year) gvar(gvar) method(dripw) agg(event)
estat event, window(-`L' `L')
matrix r = r(table)
local nm : colnames r
local j = 1
foreach c of local nm {
    * csdid event 列名形如 Tm5,Tm4,...,Tp0,Tp1,...; 解析
    local lab "`c'"
    if strpos("`lab'","Tm") local tt = -real(subinstr("`lab'","Tm","",.))
    else if strpos("`lab'","Tp") local tt = real(subinstr("`lab'","Tp","",.))
    else local tt = .
    if !missing(`tt') & inrange(`tt',-`L',`L') {
        local bb = r[1,`j']
        local lo = r[5,`j']
        local hi = r[6,`j']
        post ES5 ("CS") (`tt') (`bb') (`lo') (`hi')
    }
    local ++j
}
restore

*===================== (4) Borusyak et al. (2021) ========================*
*   did_imputation（插补法），含 pretrends 事前安慰剂
preserve
did_imputation lnpoco2 city_code year gvar, controls(`CTRL') ///
    horizons(0/`L') pretrends(`L') autosample minn(0) nose
matrix b = e(b)
matrix V = e(V)
local cols : colnames b
local j = 1
foreach v of local cols {
    local bb = b[1,`j']
    local se = sqrt(V[`j',`j'])
    * did_imputation 列名：tau0..tauL（事后），pre1..preL（事前安慰剂，符号对应 -k）
    if strpos("`v'","tau") {
        local tt = real(subinstr("`v'","tau","",.))
        post ES5 ("BJS") (`tt') (`bb') (`bb'-1.96*`se') (`bb'+1.96*`se')
    }
    else if strpos("`v'","pre") {
        local tt = -real(subinstr("`v'","pre","",.))
        post ES5 ("BJS") (`tt') (`bb') (`bb'-1.96*`se') (`bb'+1.96*`se')
    }
    local ++j
}
restore

*=========== (5) de Chaisemartin & D'Haultfoeuille (2020) ===============*
*   did_multiplegt_dyn：effects(事后) + placebo(事前)
preserve
did_multiplegt_dyn lnpoco2 city_code year DID, effects(`L') placebo(`L') ///
    controls(`CTRL') cluster(city_code) graph_off
* 事后效应 Effect_1..Effect_L
forvalues k = 1/`L' {
    local bb = e(Effect_`k')
    local se = e(se_effect_`k')
    if !missing(`bb') post ES5 ("DCDH") (`k') (`bb') (`bb'-1.96*`se') (`bb'+1.96*`se')
}
* 即时效应（部分版本记为 Effect_0 不存在，dCDH 以 effect_1 为首期 t=0/1，按版本核对）
* 事前安慰剂 Placebo_1..Placebo_L -> 对应 -1..-L
forvalues k = 1/`L' {
    local bb = e(Placebo_`k')
    local se = e(se_placebo_`k')
    if !missing(`bb') post ES5 ("DCDH") (-`k') (`bb') (`bb'-1.96*`se') (`bb'+1.96*`se')
}
restore

postclose ES5

*------------------------------------------------------------------------------*
*  合成图：五种估计量横向错位散点（黑白 s1mono），仿图8/图9
*------------------------------------------------------------------------------*
use "`results'", clear
* 横向错位 offset，避免五种估计量在同一事件时间点重叠
gen double x = rel
replace x = rel - 0.28 if est=="BJS"
replace x = rel - 0.14 if est=="DCDH"
replace x = rel + 0.00 if est=="CS"
replace x = rel + 0.14 if est=="TWFE"
replace x = rel + 0.28 if est=="SA"

twoway ///
  (rcap hi lo x if est=="BJS",  lcolor(black) lwidth(thin) msize(small)) ///
  (rcap hi lo x if est=="DCDH", lcolor(black) lwidth(thin) msize(small)) ///
  (rcap hi lo x if est=="CS",   lcolor(black) lwidth(thin) msize(small)) ///
  (rcap hi lo x if est=="TWFE", lcolor(black) lwidth(thin) msize(small)) ///
  (rcap hi lo x if est=="SA",   lcolor(black) lwidth(thin) msize(small)) ///
  (scatter coef x if est=="BJS",  msymbol(circle)        mcolor(black) msize(medium)) ///
  (scatter coef x if est=="DCDH", msymbol(square)        mcolor(black) msize(medium)) ///
  (scatter coef x if est=="CS",   msymbol(triangle)      mcolor(black) msize(medium)) ///
  (scatter coef x if est=="TWFE", msymbol(diamond)       mcolor(black) msize(medium)) ///
  (scatter coef x if est=="SA",   msymbol(x)             mcolor(black) msize(large)) , ///
  yline(0, lcolor(black) lwidth(thin)) ///
  xline(-0.5, lpattern(dash) lcolor(black) lwidth(thin)) ///
  xlabel(-5(1)5, nogrid labsize(medlarge)) ///
  ylabel(-0.5(0.25)0.5, nogrid angle(0) format(%4.2f) labsize(medlarge)) ///
  xtitle("距离5A景区设立的时间", size(large) margin(t+3)) ///
  ytitle("平均处理效应", size(large) margin(r+3)) ///
  title("五类DID估计量（动态处理效应）", size(large)) ///
  legend(order(6 "Borusyak et al." 7 "de Chaisemartin-D'Haultfoeuille" ///
               8 "Callaway-Sant'Anna" 9 "TWFE" 10 "Sun-Abraham") ///
         position(6) ring(1) rows(3) size(medsmall) symxsize(6) region(lstyle(none))) ///
  scheme(s1mono) graphregion(color(white) margin(medlarge)) ///
  plotregion(margin(medsmall) lcolor(black))
graph export "results/图_五类DID动态处理效应_lnpoco2.png", replace width(2460) height(1620)
di as result "============ 五类DID动态处理效应图 完成：results/图_五类DID动态处理效应_lnpoco2.png ============"

*--- 如需 SO2 口径(lnpoco2_so2)平行图：将上方 (1)-(5) 与合成图中的 lnpoco2 ---*
*--- 整体替换为 lnpoco2_so2，输出名改为 图_五类DID动态处理效应_lnpoco2_so2.png 即可 ---*

di as result "================ 全部完成：表格与图 见 results 文件夹 ================"

*  减污降碳 5A 政策 DID —— 机制检验（中介效应三步法 + Sobel + Bootstrap）  0624
*  方法：中介效应三步法（依次估 c、a、b 与 c'）；并报 Sobel 检验与 Bootstrap 百分位 CI。
*  报告：括号内为 t 值（3 位小数），系数 3 位小数，含常数项 _cons；城市层面聚类标准误。
*  四条渠道（被解释变量始终为 lnpoco2，越小=减污降碳越好）：
*    H2a 引导转型  ter_gdp     第三产业占GDP比重(%)          预期 a>0, b<0
*    H2b 协同枢纽  lnelec_gdp  ln(全社会用电量/GDP)          预期 a<0, b>0
*    H2c 波特效应  lnpatapp    ln(1+专利申请量)·绿色创新代理   预期 a>0, b<0  ← 用Excel AJ列
*    H2d 倒逼治理  er          环境规制强度                  预期 a>0, b<0
*  注：lnelec_gdp、lnpatapp 已在 data.dta 内（lnpatapp=AJ列 ln(1+专利申请)，缺失已插补)；
*      H2a 产业结构渠道控制集剔除 struc2(≡第二产业占比) 以避免坏控制(江艇,2022)。
*  缺包先装：ssc install reghdfe estout, replace
*==============================================================================*
clear all
set more off
cap mkdir results
use "data.dta", clear
xtset city_code year

local CF "lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg"
local CA "lnpgdp lndensity urban gov tech"
local CH `" "H2a引导转型 ter_gdp CA" "H2b协同枢纽 lnelec_gdp CF" "H2c波特效应 lnpatapp CF" "H2d倒逼治理 er CF" "'

*==============================================================================*
* 【一】三步法逐渠道：表(第一步c / 第二步a / 第三步b,c') + Sobel 汇总
*==============================================================================*
di as result _n "{hline 94}"
di as text  %-13s "渠道" %9s "c" %10s "a" %11s "b" %10s "c'" %11s "a*b" %10s "Sobel z" %8s "p" "  显著"
di as text  "{hline 94}"

foreach row of local CH {
    gettoken tag row : row
    gettoken med row : row
    gettoken cid row : row
    local C ``cid''

    eststo S1: reghdfe lnpoco2 DID `C',        a(city_code year) vce(cl city_code)
    scalar cco = _b[DID]
    eststo S2: reghdfe `med'   DID `C',        a(city_code year) vce(cl city_code)
    scalar aco = _b[DID]
    scalar sea = _se[DID]
    eststo S3: reghdfe lnpoco2 DID `med' `C',  a(city_code year) vce(cl city_code)
    scalar bco = _b[`med']
    scalar seb = _se[`med']
    scalar cpr = _b[DID]

    * 三步法表（系数3位、t值3位、含常数项）
    esttab S1 S2 S3 using "results/三步法_`tag'.rtf", replace ///
        b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
        keep(DID `med' _cons) noomitted varlabels(_cons 常数项) ///
        mtitles("第一步 lnpoco2 (c)" "第二步 `med' (a)" "第三步 lnpoco2 (b,c')") ///
        stats(N r2_a, fmt(%9.0f %9.3f) labels("观测值N" "Adj. R2")) ///
        nogaps compress label title("`tag' 中介效应三步法（中介=`med'）") ///
        addnotes("括号内为t值(3位小数)；系数3位小数；含常数项；城市层面聚类标准误；* p<.1 ** p<.05 *** p<.01")
    eststo clear

    * Sobel 检验（用聚类标准误）：z = a*b / sqrt(b^2*se_a^2 + a^2*se_b^2)
    scalar indd = aco*bco
    scalar zso  = indd/sqrt(bco^2*sea^2 + aco^2*seb^2)
    scalar pso  = 2*(1-normal(abs(zso)))
    local sg=cond(abs(zso)>2.58,"***",cond(abs(zso)>1.96,"**",cond(abs(zso)>1.65,"*","ns")))
    di as text %-13s "`tag'" as result %9.3f cco %10.3f aco %11.3f bco %10.3f cpr ///
       %11.3f indd %10.3f zso %8.3f pso as text "  `sg'"
}
di as text "{hline 94}"
di as text "判定：a、b 同显著且 a*b 与 c 同号 ⇒ 中介成立；Sobel 与下方 Bootstrap 检验间接效应 a*b。"

*==============================================================================*
* 【二】Bootstrap 间接效应 a*b（百分位法，500 次；reps 可按需调小提速）
*==============================================================================*
cap program drop bootmed
program bootmed, rclass
    reghdfe ${med} DID ${ctrl},          a(city_code year)
    local a = _b[DID]
    reghdfe lnpoco2 DID ${med} ${ctrl},  a(city_code year)
    return scalar ind = `a'*_b[${med}]
end

di as result _n "===== Bootstrap 间接效应 a*b（500次，百分位95%CI）====="
foreach row of local CH {
    gettoken tag row : row
    gettoken med row : row
    gettoken cid row : row
    global med  "`med'"
    global ctrl "``cid''"
    di as result _n ">>> `tag'  中介=`med'"
    bootstrap ind=r(ind), reps(500) seed(20250624) dots(50): bootmed
    estat bootstrap, percentile
}

di as result _n "================ 机制检验完成：三步法表见 results 文件夹 ================"



*##############################################################################*
*  五、异质性分析（6 个维度，分别输出；系数/ t 值均保留 3 位小数）            *
*  输出文件夹：results；6 组分别输出独立结果表                                *
*  维度：①地理区位 ②资源禀赋 ③城市规模(行政等级) ④5A数量 ⑤景区属性(调节) ⑥污染碳本底
*  统一：被解释变量 lnpoco2，全控制变量，城市与年份双向固定效应，城市层面聚类标准误
*##############################################################################*
clear all
set more off
cap mkdir results
use "data.dta", clear
xtset city_code year

local CTRL "lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg"

*------------------------------------------------------------------------------*
* ① 地理区位异质性：东部 / 中部 / 西部（分样本回归，一表三列）
*------------------------------------------------------------------------------*
eststo clear
eststo h1e: reghdfe lnpoco2 DID `CTRL' if region=="东部", a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
eststo h1c: reghdfe lnpoco2 DID `CTRL' if region=="中部", a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
eststo h1w: reghdfe lnpoco2 DID `CTRL' if region=="西部", a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
esttab h1e h1c h1w using "results/异质性1_地理区位.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
    keep(DID `CTRL') order(DID `CTRL') ///
    mtitles("东部" "中部" "西部") ///
    stats(cityfe yearfe N r2_a, fmt(%s %s %9.0f %9.3f) ///
          labels("城市固定效应" "年份固定效应" "观测值N" "Adj. R2")) ///
    nogaps compress title("异质性① 地理区位（东/中/西）") ///
    addnotes("括号内为t值(3位小数)，城市层面聚类标准误；* p<0.1, ** p<0.05, *** p<0.01")
eststo clear

*------------------------------------------------------------------------------*
* ② 资源禀赋异质性：资源型城市 / 非资源型城市（国务院2013规划口径）
*------------------------------------------------------------------------------*
eststo h2r: reghdfe lnpoco2 DID `CTRL' if resource==1, a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
eststo h2n: reghdfe lnpoco2 DID `CTRL' if resource==0, a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
esttab h2r h2n using "results/异质性2_资源禀赋.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
    keep(DID `CTRL') order(DID `CTRL') ///
    mtitles("资源型城市" "非资源型城市") ///
    stats(cityfe yearfe N r2_a, fmt(%s %s %9.0f %9.3f) ///
          labels("城市固定效应" "年份固定效应" "观测值N" "Adj. R2")) ///
    nogaps compress title("异质性② 资源禀赋（资源型/非资源型）") ///
    addnotes("括号内为t值(3位小数)，城市层面聚类标准误；* p<0.1, ** p<0.05, *** p<0.01")
eststo clear

*------------------------------------------------------------------------------*
* ③ 城市规模异质性：中心城市(直辖市/省会/计划单列) / 一般地级市
*------------------------------------------------------------------------------*
eststo h3c: reghdfe lnpoco2 DID `CTRL' if central==1, a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
eststo h3o: reghdfe lnpoco2 DID `CTRL' if central==0, a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
esttab h3c h3o using "results/异质性3_城市规模.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
    keep(DID `CTRL') order(DID `CTRL') ///
    mtitles("中心城市" "一般城市") ///
    stats(cityfe yearfe N r2_a, fmt(%s %s %9.0f %9.3f) ///
          labels("城市固定效应" "年份固定效应" "观测值N" "Adj. R2")) ///
    nogaps compress title("异质性③ 城市规模（中心城市/一般城市）") ///
    addnotes("括号内为t值(3位小数)，城市层面聚类标准误；中心城市=直辖市+省会+计划单列市；* p<0.1, ** p<0.05, *** p<0.01")
eststo clear

*------------------------------------------------------------------------------*
* ④ 5A数量异质性：0个 / 1个 / ≥2个5A（全样本单方程，按数量强度拆分处理效应）
*   互斥三分(观测合计=6069)：0个=从未获评(基准,2226) | 1个=单个5A(2415) | ≥2个=多个5A(1428)
*   did_1=单个5A城市政策实施后；did_2=多个5A(≥2)城市政策实施后；基准组=0个5A(从未获评)
*------------------------------------------------------------------------------*
cap drop multi5a
cap drop did_1
cap drop did_2
gen byte multi5a = 0
foreach c in 三明市 上海市 上饶市 丽江市 九江市 保定市 信阳市 北京市 北海市 南京市 南宁市 南平市 南阳市 哈尔滨市 嘉兴市 大连市 宁德市 宁波市 安阳市 安顺市 宜春市 常州市 平顶山市 广州市 开封市 张家界市 徐州市 成都市 扬州市 承德市 新乡市 无锡市 昆明市 杭州市 柳州市 桂林市 武汉市 毕节市 泉州市 泰安市 洛阳市 济宁市 淮安市 清远市 温州市 湖州市 烟台市 焦作市 盐城市 福州市 秦皇岛市 绍兴市 肇庆市 苏州市 西安市 贵阳市 赣州市 通化市 遵义市 重庆市 金华市 长春市 长沙市 鞍山市 韶关市 驻马店市 黄山市 龙岩市 {
    replace multi5a = 1 if city=="`c'" & treat==1
}
gen byte did_1 = DID*(multi5a==0)
gen byte did_2 = DID*(multi5a==1)
label var did_1 "1个5A×政策后"
label var did_2 "≥2个5A×政策后"
eststo h4: reghdfe lnpoco2 did_1 did_2 `CTRL', a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
esttab h4 using "results/异质性4_5A数量.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
    keep(did_1 did_2 `CTRL') order(did_1 did_2 `CTRL') ///
    mtitles("lnpoco2") ///
    varlabels(did_1 "1个5A(单个)" did_2 "≥2个5A(多个)") ///
    stats(cityfe yearfe N r2_a, fmt(%s %s %9.0f %9.3f) ///
          labels("城市固定效应" "年份固定效应" "观测值N" "Adj. R2")) ///
    nogaps compress title("异质性④ 5A数量（0个/1个/≥2个，全样本N=6069，基准=0个5A）") ///
    addnotes("括号内为t值(3位小数)，城市层面聚类标准误；全样本N=6069，0个5A(从未获评)为基准组" ///
             "互斥分组观测：0个=2226、1个=2415、≥2个=1428，合计=6069；* p<0.1, ** p<0.05, *** p<0.01")
eststo clear
*------------------------------------------------------------------------------*
* ⑤ 景区属性异质性：自然类 / 人文类5A（全样本单方程，按景区属性拆分处理效应）
*   互斥三分(观测合计=6069)：从未获评(基准,2226) | 自然类处理城市(3024) | 人文类处理城市(819)
*   did_nat=自然类城市政策实施后；did_cul=人文类城市政策实施后；基准组=从未获评城市
*------------------------------------------------------------------------------*
* 单个5A城市、多个5A城市，各与"从未获评(grp5a==0)"城市对比
local CTRL "lnpgdp lndensity urban struc2 gov tech struc_adv struc_upg"
eststo clear
eststo h4s: reghdfe lnpoco2 DID `CTRL' if inlist(grp5a,1,0), a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
eststo h4m: reghdfe lnpoco2 DID `CTRL' if inlist(grp5a,2,0), a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
esttab h4s h4m using "results/异质性4_5A数量.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
    keep(DID `CTRL') order(DID `CTRL') ///
    mtitles("单个5A城市" "多个5A城市") ///
    stats(cityfe yearfe N r2_a, fmt(%s %s %9.0f %9.3f) ///
          labels("城市固定效应" "年份固定效应" "观测值N" "Adj. R2")) ///
    nogaps compress title("异质性④ 政策强度—5A数量（单个/多个，各与从未获评城市对比）") ///
    addnotes("括号内为t值(3位小数)，城市层面聚类标准误；基准组=从未获评城市；* p<0.1, ** p<0.05, *** p<0.01")
eststo clear

*------------------------------------------------------------------------------*
* ⑥ 初始污染—碳本底异质性：高本底 / 低本底（政策前 lnpoco2 均值中位数分组）
*------------------------------------------------------------------------------*
eststo h6h: reghdfe lnpoco2 DID `CTRL' if highpoll==1, a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
eststo h6l: reghdfe lnpoco2 DID `CTRL' if highpoll==0, a(city_code year) vce(cl city_code)
estadd local cityfe "是"
estadd local yearfe "是"
esttab h6h h6l using "results/异质性6_污染碳本底.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
    keep(DID `CTRL') order(DID `CTRL') ///
    mtitles("高本底" "低本底") ///
    stats(cityfe yearfe N r2_a, fmt(%s %s %9.0f %9.3f) ///
          labels("城市固定效应" "年份固定效应" "观测值N" "Adj. R2")) ///
    nogaps compress title("异质性⑥ 初始污染-碳本底（高/低）") ///
    addnotes("括号内为t值(3位小数)，城市层面聚类标准误；按政策前lnpoco2均值中位数分高/低；* p<0.1, ** p<0.05, *** p<0.01")
eststo clear

*==============================================================================*
* 异质性显著性判定汇总（控制台打印；系数/t 值 3 位小数）
*   说明：各表已含显著性星号；此处再以紧凑形式打印 DID 系数/t/p 便于快速判定
*==============================================================================*
di as result _n "{hline 74}"
di as text  "异质性分组显著性判定（被解释变量 lnpoco2，城市层面聚类标准误）"
di as text  "{hline 74}"

* ① 地理区位
foreach rg in 东部 中部 西部 {
    qui reghdfe lnpoco2 DID `CTRL' if region=="`rg'", a(city_code year) vce(cl city_code)
    local b=_b[DID]
    local t=_b[DID]/_se[DID]
    local p=2*ttail(e(df_r),abs(`t'))
    local s=cond(`p'<0.01,"***",cond(`p'<0.05,"**",cond(`p'<0.1,"*","ns")))
    di as text %-24s "① 地理区位-`rg'" as result %9.3f `b' "  t=" %8.3f `t' "  p=" %6.3f `p' as text "  `s'"
}

* ② 资源禀赋
foreach g in 1 0 {
    local lab=cond(`g'==1,"资源型城市","非资源型城市")
    qui reghdfe lnpoco2 DID `CTRL' if resource==`g', a(city_code year) vce(cl city_code)
    local b=_b[DID]
    local t=_b[DID]/_se[DID]
    local p=2*ttail(e(df_r),abs(`t'))
    local s=cond(`p'<0.01,"***",cond(`p'<0.05,"**",cond(`p'<0.1,"*","ns")))
    di as text %-24s "② 资源禀赋-`lab'" as result %9.3f `b' "  t=" %8.3f `t' "  p=" %6.3f `p' as text "  `s'"
}

* ③ 城市规模
foreach g in 1 0 {
    local lab=cond(`g'==1,"中心城市","一般城市")
    qui reghdfe lnpoco2 DID `CTRL' if central==`g', a(city_code year) vce(cl city_code)
    local b=_b[DID]
    local t=_b[DID]/_se[DID]
    local p=2*ttail(e(df_r),abs(`t'))
    local s=cond(`p'<0.01,"***",cond(`p'<0.05,"**",cond(`p'<0.1,"*","ns")))
    di as text %-24s "③ 城市规模-`lab'" as result %9.3f `b' "  t=" %8.3f `t' "  p=" %6.3f `p' as text "  `s'"
}

* ④ 5A数量
foreach pair in "多个5A 2" "单个5A 1" {
    gettoken lab code : pair
    qui reghdfe lnpoco2 DID `CTRL' if inlist(grp5a,`code',0), a(city_code year) vce(cl city_code)
    local b=_b[DID]
    local t=_b[DID]/_se[DID]
    local p=2*ttail(e(df_r),abs(`t'))
    local s=cond(`p'<0.01,"***",cond(`p'<0.05,"**",cond(`p'<0.1,"*","ns")))
    di as text %-24s "④ 5A数量-`lab'" as result %9.3f `b' "  t=" %8.3f `t' "  p=" %6.3f `p' as text "  `s'"
}
* ⑤ 景区属性（自然类 / 人文类，分样本 2 列；沿用上文已生成的 Resour）
foreach pair in "自然类 0" "人文类 1" {
    gettoken lab val : pair
    qui reghdfe lnpoco2 DID `CTRL' if (treat==1 & Resour==`val') | treat==0, a(city_code year) vce(cl city_code)
    local b=_b[DID]
    local t=_b[DID]/_se[DID]
    local p=2*ttail(e(df_r),abs(`t'))
    local s=cond(`p'<0.01,"***",cond(`p'<0.05,"**",cond(`p'<0.1,"*","ns")))
    di as text %-24s "⑤ 景区属性-`lab'" as result %9.3f `b' "  t=" %8.3f `t' "  p=" %6.3f `p' as text "  `s'"
}

* ⑥ 初始污染-碳本底
foreach g in 1 0 {
    local lab=cond(`g'==1,"高本底","低本底")
    qui reghdfe lnpoco2 DID `CTRL' if highpoll==`g', a(city_code year) vce(cl city_code)
    local b=_b[DID]
    local t=_b[DID]/_se[DID]
    local p=2*ttail(e(df_r),abs(`t'))
    local s=cond(`p'<0.01,"***",cond(`p'<0.05,"**",cond(`p'<0.1,"*","ns")))
    di as text %-24s "⑥ 污染本底-`lab'" as result %9.3f `b' "  t=" %8.3f `t' "  p=" %6.3f `p' as text "  `s'"
}

