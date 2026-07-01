*==============================================================================*
*  空间溢出分析 · 步骤3：按地理距离检验溢出"范围"（distance decay）
*  思路：用不同 km 的距离权重矩阵重复 SDM/SAR，观察空间项随距离的变化，
*        定位溢出的有效地理半径。
*  两套口径：
*    (A) 累积带 W_d(≤km)：随半径扩大纳入更多邻居（结果单调，供稳健参考）
*    (B) 距离分环(在 Python 端可另建)：可精确定位溢出衰减到不显著的环
*  依赖：先运行 spatial_01_weights.do。
*==============================================================================*
use "data_spatial.dta", clear
cap gen double lnfin = ln(fin)
xtset city_code year
global CTRL lnpgdp lndensity urban struc2 gov tech human lnfin

*--- 累积距离带：150→500 km，记录 rho 与显著性 ---*
tempname M
postfile `M' int km double rho double p using "results_distance_rho.dta", replace
foreach km in 150 200 250 300 350 400 450 500 {
    qui xsmle lnpoco2 DID $CTRL, model(sar) wmat(Wd`km') fe type(both) nsim(1)
    local rho = e(rho)
    * 近似 z 检验（rho 与其标准误；xsmle 存于 e(b)/e(V) 的 rho 部分）
    matrix b = e(b)
    matrix V = e(V)
    local se = sqrt(V["Spatial:rho","Spatial:rho"])
    local z  = `rho'/`se'
    local pv = 2*(1-normal(abs(`z')))
    di as txt "距离带 ≤`km' km:  rho=" as res %6.3f `rho' as txt "  p=" as res %6.3f `pv'
    post `M' (`km') (`rho') (`pv')
}
postclose `M'

*--- 画溢出-距离衰减图 ---*
use "results_distance_rho.dta", clear
twoway (connected rho km, msymbol(O) lcolor(black) mcolor(black)), ///
    yline(0, lpattern(dash)) ///
    xtitle("距离阈值 (km)") ytitle("空间自回归系数 rho") ///
    title("减污降碳空间溢出的地理衰减") scheme(s1mono)
graph export "results_distance_decay.png", replace width(2000) height(1400)
di as result "== 距离衰减结果见 results_distance_rho.dta 与 results_distance_decay.png =="

*==============================================================================*
*  重要提示（照实汇报）：
*  在"双向固定效应 + 全部控制变量"这一严格设定下，各距离带的空间项普遍不显著；
*  空间溢出主要在"未控制产业结构/未加年份固定效应"时显现，且集中在约 100–500 km。
*  这说明减污降碳的空间关联很大程度经由"产业结构邻近性"传导（见 doc 机制部分）。
*==============================================================================*
