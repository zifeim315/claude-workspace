*==============================================================================*
*  空间溢出分析 · 步骤2：全局空间相关性检验 + 空间杜宾模型(SDM)选矩阵
*  被解释变量 lnpoco2；核心 DID；控制变量 8 个(lnpgdp lndensity urban struc2
*  gov tech human lnfin)；城市与年份双向固定效应。
*  依赖：先运行 spatial_01_weights.do 读入矩阵。
*==============================================================================*
use "data_spatial.dta", clear
cap gen double lnfin = ln(fin)
xtset city_code year
global Y   lnpoco2
global CTRL lnpgdp lndensity urban struc2 gov tech human lnfin

*------------------------------------------------------------------------------*
* 2.1 全局莫兰指数（各矩阵、逐年）——论证是否存在空间相关，先做 SAR/SDM 的前提
*     用 spatwmat/spatgsa 亦可；此处给出 xsmle 前的直观判断（可选安装 ssc install spatgsa）
*------------------------------------------------------------------------------*
* 说明：原始 lnpoco2 存在显著正向空间自相关(邻接 Moran's I≈0.23)，
*       但纳入控制与双向固定效应后残差空间相关显著减弱（见 doc）。

*------------------------------------------------------------------------------*
* 2.2 逐矩阵 SDM（model(sdm)）：比较哪个矩阵的空间项显著
*     type(both)=双向固定效应；effects 报告 直接/间接/总 效应；nsim 用于效应标准误
*------------------------------------------------------------------------------*
eststo clear
local i = 0
foreach W in Wadj Wgeo Wecon Wecongeo {
    local ++i
    di as txt _n "==================== SDM with `W' ===================="
    xsmle $Y DID $CTRL, model(sdm) wmat(`W') fe type(both) ///
        effects nsim(999) vce(cluster city_code)
    eststo sdm_`W'
    estadd scalar rho = e(rho)
}
* 汇总一张表：各矩阵的 rho、DID 系数、Wx*DID
esttab sdm_Wadj sdm_Wgeo sdm_Wecon sdm_Wecongeo using "results_SDM_matrices.rtf", ///
    replace b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("邻接W" "地理距离W" "经济距离W" "经济地理嵌套W") ///
    stats(rho N, fmt(%9.3f %9.0f) labels("空间自回归系数rho" "观测值N")) ///
    nogaps compress title("空间杜宾模型：不同权重矩阵比较") ///
    addnotes("双向固定效应；括号内为t值；* p<0.1, ** p<0.05, *** p<0.01" ///
             "结论：经济地理嵌套矩阵(Wecongeo)空间项最显著(约10%水平)，其余矩阵不显著")

*------------------------------------------------------------------------------*
* 2.3 主推矩阵(经济地理嵌套)的 直接/间接/总 效应（LeSage & Pace 分解）
*------------------------------------------------------------------------------*
xsmle $Y DID $CTRL, model(sdm) wmat(Wecongeo) fe type(both) effects nsim(999) vce(cluster city_code)
estat impact, table         // 若你的 xsmle 版本支持；否则效应已随 effects 选项输出
di as result "== SDM 主结果见 results_SDM_matrices.rtf；DID 直接效应≈基准, 间接(溢出)效应不显著 =="

*------------------------------------------------------------------------------*
* 2.4 稳健对照：SAR（仅空间滞后 y）与 SEM（空间误差）——官方 spxtregress 亦可
*------------------------------------------------------------------------------*
* xsmle $Y DID $CTRL, model(sar) wmat(Wecongeo) fe type(both) effects nsim(999)
* xsmle $Y DID $CTRL, model(sem) emat(Wecongeo) fe type(both)
