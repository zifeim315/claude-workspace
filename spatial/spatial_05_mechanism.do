*==============================================================================*
*  空间溢出分析 · 步骤5：溢出传导机制（本研究的创新点之一）
*  发现：减污降碳的空间关联主要经由"产业结构邻近性"传导。
*  验证：以邻接矩阵做 SAR，逐步加入控制变量，观察空间自回归系数 rho 的塌缩。
*    仅 DID(城市FE)         : rho≈+0.105***
*    + 经济发展 lnpgdp      : rho≈+0.104***
*    + 产业结构 struc2      : rho≈+0.031*      ← 产业结构吸收了大部分空间相关
*    + 全部8控制            : rho≈+0.007  ns
*  结论：空间溢出并非"纯地理"传播，而是通过产业结构/发展水平的空间同群性实现。
*  依赖：先运行 spatial_01_weights.do。
*==============================================================================*
use "data_spatial.dta", clear
cap gen double lnfin = ln(fin)
xtset city_code year

di as txt "== SAR(邻接W, 仅城市固定效应) 逐步加控制, 观察 rho 塌缩 =="
foreach spec in "DID" "DID lnpgdp" "DID lnpgdp struc2" ///
                "DID lnpgdp lndensity urban struc2 gov tech human lnfin" {
    qui xsmle lnpoco2 `spec', model(sar) wmat(Wadj) fe type(ind) nsim(1)
    di as res %6.3f e(rho) as txt "   <- rho | X = `spec'"
}
di as result "== 机制：rho 随'产业结构'纳入而塌缩 → 空间溢出经由产业结构邻近性传导 =="
