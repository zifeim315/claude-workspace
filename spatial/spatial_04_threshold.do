*==============================================================================*
*  空间溢出分析 · 步骤4：门槛 / 分段异质性溢出（SDM + 门槛）
*  两种做法：
*   (A) 面板门槛模型(Hansen 1999, xthreg)：以 经济发展水平/产业结构 为门槛变量，
*       检验"本地效应 + 空间溢出强度"是否存在结构性分段。
*   (B) 分区制 SDM：按门槛把样本分高/低两组，或在 SDM 中让 W×DID 随门槛虚拟变量分段，
*       比较不同区制下溢出系数差异（分段异质性溢出）。
*  依赖：先运行 spatial_01_weights.do。需安装：ssc install xthreg
*==============================================================================*
use "data_spatial.dta", clear
cap gen double lnfin = ln(fin)
xtset city_code year
global CTRL lnpgdp lndensity urban struc2 gov tech human lnfin

*------------------------------------------------------------------------------*
* (A) Hansen 面板门槛：门槛变量=经济发展水平 lnpgdp（也可换 struc2 产业结构）
*     被解释 lnpoco2，核心 DID，检验单一/双重门槛
*------------------------------------------------------------------------------*
xthreg lnpoco2 DID $CTRL, rx(DID) qx(lnpgdp) thnum(2) grid(300) trim(0.01 0.01) bs(300 300)
* 解读：若门槛显著，说明政策的减污降碳效应在不同发展水平区间存在分段差异。
estimates store thr_dev

*------------------------------------------------------------------------------*
* (B) 分区制 SDM：按 lnpgdp 城市均值中位数分高/低发展两组，
*     在经济地理嵌套矩阵下让政策空间溢出项(W×DID)分区制估计
*------------------------------------------------------------------------------*
* 生成城市均值发展分组
bysort city_code: egen dev_m = mean(lnpgdp)
qui sum dev_m, detail
gen byte hidev = dev_m > r(p50)
* 构造空间滞后 W×DID（用 Mata 从 Wecongeo 生成），再按区制交互
mata:
    W = st_matrix("Wecongeo")
    st_view(y=., ., "DID")
    id = st_data(., "city_code"); yr = st_data(., "year")
    // 逐年做 W*DID
    info = panelsetup(yr, 1)
end
* 简便起见：在 xsmle 的 SDM 中直接纳入分组交互(近似分段溢出)
gen DID_hi = DID*hidev
gen DID_lo = DID*(1-hidev)
xsmle lnpoco2 DID_lo DID_hi $CTRL, model(sdm) wmat(Wecongeo) fe type(both) ///
    effects nsim(999) vce(cluster city_code)
estimates store sdm_regime
di as result "== 门槛/分区制结果：thr_dev(Hansen门槛) 与 sdm_regime(分段溢出) =="

*==============================================================================*
*  提示：本研究数据在双向固定效应下政策空间溢出整体不显著，
*  门槛/分段主要用于刻画"本地减排效应"的区制差异（高/低发展、产业结构门槛）。
*==============================================================================*
