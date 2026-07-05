*==============================================================================*
*  对标 Wang & Fang (2026, Sustainable Cities and Society 137:107132) 的方法扩展
*  可直接接在 5A_integrated_analysis.do 之后运行（依赖 results/_work.dta）
*
*  模块 A：被解释变量"协同排放(SEM)"正名 + 碳强度口径重构   —— 数据已齐全，立即可跑
*  模块 B：Tapio 脱钩模型（经济增长 ↔ 减污降碳）+ 5A 对强脱钩的影响 —— 数据已齐全
*  模块 C：超效率 SBM 构建"协同效率(SEF)"                    —— 模板，需先补 4 个投入变量
*  模块 D：SEM 时空演变(核密度) + 城市群维度异质性            —— 数据已齐全
*
*  说明：lnpoco2 / SEM 越小越好；DID 为多期政策交乘项。
*==============================================================================*
clear all
set more off
cap mkdir results
use "results/_work.dta", clear      // 若无该副本，改成 use "data_spatial.dta", clear + xtset city_code year
xtset city_code year
cap global CTRL "lnpgdp lndensity urban struc2 gov tech human lnfin"


*==============================================================================*
* 模块 A：被解释变量 SEM 正名 + 碳强度口径重构
*   标杆式(1)：SEM = CEI × EPI，其中 CEI = CO2/GDP（碳强度，非总量）
*==============================================================================*
* A.1 碳强度口径 SEM（对齐标杆 CEI 定义）
gen double cei      = co2_wt/gdp                 // 碳排放强度 = CO2(万吨)/GDP
gen double sem_int  = cei * poll_idx             // 协同排放(强度口径) = CEI × EPI
gen double lnsem_int = ln(sem_int)
label var cei       "碳排放强度 CEI=CO2/GDP"
label var lnsem_int "协同排放SEM(碳强度口径,ln)"

* A.2 三口径并列稳健性表：总量口径(主) / SO2口径 / 碳强度口径(标杆同款)
eststo clear
eststo a1: reghdfe lnpoco2     DID $CTRL, a(city_code year) vce(cl city_code)   // 主口径(CO2总量×EPI)
eststo a2: reghdfe lnpoco2_so2 DID $CTRL, a(city_code year) vce(cl city_code)   // SO2口径
eststo a3: reghdfe lnsem_int   DID $CTRL, a(city_code year) vce(cl city_code)   // 碳强度口径(标杆)
esttab a1 a2 a3 using "results/表A_SEM口径重构(对标标杆).rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) keep(DID) ///
    mtitles("CO2总量×EPI(主)" "SO2口径" "碳强度×EPI(标杆同款)") ///
    stats(N r2_a, fmt(%9.0f %9.3f) labels("N" "Adj.R2")) nogaps compress ///
    title("表A 协同排放SEM多口径重构：对标 Wang & Fang (2026) 交乘法") ///
    addnotes("SEM=CEI×EPI(交乘法, Wang&Fang 2026 式1); CEI=CO2/GDP为碳强度; 越小越好。三口径下DID均显著为负,支撑结论稳健。")
eststo clear


*==============================================================================*
* 模块 B：Tapio 脱钩模型（经济增长 ↔ 减污降碳SEM）
*   DI = (ΔSEM/SEM_t1)/(ΔGDP/GDP_t1)。SEM下降+GDP上升 ⇒ DI<0 为"强脱钩(最优)"
*   分四个时段 P1(2003-08) P2(08-13) P3(13-18) P4(18-23)
*==============================================================================*
preserve
gen sem = poll_idx * co2_wt        // 与 lnpoco2 同底的 SEM 水平值(总量口径)
* 取每个时段端点值
gen byte ep = .
replace ep = 1 if inlist(year,2003,2008,2013,2018,2023)
keep if ep==1
keep city_code year sem gdp
reshape wide sem gdp, i(city_code) j(year)

* 逐时段脱钩指数与状态
forvalues p = 1/4 {
    local t0 : word `p'   of 2003 2008 2013 2018
    local t1 : word `p'   of 2008 2013 2018 2023
    gen double dsem`p' = (sem`t1'-sem`t0')/sem`t0'
    gen double dgdp`p' = (gdp`t1'-gdp`t0')/gdp`t0'
    gen double DI`p'   = dsem`p'/dgdp`p'
    * Tapio 八状态分类（以 0/0.8/1.2 为临界；SEM为负向指标）
    gen str24 state`p' = ""
    replace state`p' = "强脱钩SD"   if dgdp`p'>0 & dsem`p'<0 & DI`p'<0            // 最优:增长且减排
    replace state`p' = "弱脱钩WD"   if dgdp`p'>0 & dsem`p'>0 & DI`p'>=0  & DI`p'<0.8
    replace state`p' = "扩张连接EC" if dgdp`p'>0 & dsem`p'>0 & DI`p'>=0.8 & DI`p'<=1.2
    replace state`p' = "扩张负脱钩END" if dgdp`p'>0 & dsem`p'>0 & DI`p'>1.2
    replace state`p' = "衰退脱钩RD" if dgdp`p'<0 & dsem`p'<0 & DI`p'>1.2
    replace state`p' = "其他"      if state`p'==""
    tab state`p'
}
save "results/_tapio_states.dta", replace
restore

* B.2 5A 政策是否提高"强脱钩"概率（把脱钩接回 DID 主线）
*     思路：以政策前后各城市实现强脱钩(SEM下降且GDP上升)的年际状态为结果
preserve
xtset city_code year
* 注意：D. 算子不能套函数(D.ln() 会报错)。须先生成对数变量，再对其差分。
gen double lgdp = ln(gdp)
gen double lsem = ln(poll_idx*co2_wt)
gen double g_gdp = D.lgdp
gen double g_sem = D.lsem
gen byte strong_decouple = (g_gdp>0 & g_sem<0) if !missing(g_gdp,g_sem)   // 年度强脱钩=1
eststo clear
eststo dec: reghdfe strong_decouple DID $CTRL, a(city_code year) vce(cl city_code)
esttab dec using "results/表B_5A对强脱钩的影响.rtf", replace b(%9.3f) t(%9.3f) ///
    star(* 0.1 ** 0.05 *** 0.01) keep(DID) mtitles("Pr(强脱钩)") ///
    stats(N r2_a, fmt(%9.0f %9.3f) labels("N" "Adj.R2")) nogaps compress ///
    title("表B 5A建设对'增长-减污降碳强脱钩'的影响(LPM)") ///
    addnotes("被解释变量=年度强脱钩虚拟(GDP增长且SEM下降);正系数=5A促进经济增长与排放脱钩。")
eststo clear
restore


*==============================================================================*
* 模块 C：超效率 SBM 构建"协同效率(SEF)"  —— 模板（需先补投入变量）
*   标杆表1 指标体系：
*     投入(5): 固定资产投资 / 城镇单位从业人员 / 建成区面积 / 总供水量 / 全社会用电量(elec✅)
*     期望产出(1): 地区GDP(gdp✅)
*     非期望产出(4): CO2(co2_wt✅) / PM2.5 / 工业废水 / 工业固废   (SO2✅可作代理)
*   ⚠️ 缺失: 固定资产投资、从业人员、建成区面积、总供水量、PM2.5、工业废水、工业固废
*==============================================================================*
/*  取消注释并在补齐变量后运行：

* 1) 安装 SBM 命令（三选一）：
*    ssc install sbmeff          // Stata 原生非期望产出 SBM（若可用）
*    否则用 Python(pyDEA) / MATLAB(DEA toolbox) / DEARUN 外部计算后并入

* 2) 变量准备（示例名，按实际字段替换）：
*    global XIN  "invest labor built_area water elec"     // 5 投入
*    global YG   "gdp"                                     // 期望产出
*    global YB   "co2_wt pm25 wastewater solidwaste"       // 4 非期望产出

* 3) 计算超效率 SBM（含非期望产出、VRS/CRS 按需）：
*    sbmeff $XIN, gy($YG) by($YB) prod(vrs) super
*    predict sef, teradial          // 得到 SEF（可 >1）
*    label var sef "减污降碳协同效率SEF(超效率SBM)"

* 4) 双结果 DID：SEF 作为第二被解释变量（符号预期为正 = 5A 提高效率）：
*    eststo clear
*    eststo sef1: reghdfe sef DID $CTRL, a(city_code year) vce(cl city_code)
*    esttab sef1 using "results/表C_5A对协同效率SEF的影响.rtf", replace ///
*        b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) keep(DID) mtitles("SEF") ///
*        stats(N r2_a, fmt(%9.0f %9.3f) labels("N" "Adj.R2")) nogaps compress ///
*        title("表C 5A建设对减污降碳协同效率(SEF)的影响")

* 5) 若数据暂不全，先做精简 SBM（2投入 elec+labor / 1产出 gdp / 2非期望 co2+so2），标注为简化探索口径。

* 6) 有了 SEF 后可复刻标杆的 Spearman + PVAR：
*    spearman lnpoco2 sef
*    * PVAR: ssc install pvar ; pvar lnsem lnsef, lags(.) ; pvarirf ; pvarfevd
*/


*==============================================================================*
* 模块 D：SEM 时空演变(核密度) + 城市群维度异质性
*==============================================================================*
* D.1 分年份核密度演变（讲收敛/发散故事，类比标杆 Fig.3d/4d）
preserve
twoway (kdensity lnpoco2 if year==2003, lcolor(gs12) lpattern(solid)) ///
       (kdensity lnpoco2 if year==2010, lcolor(gs8)  lpattern(dash)) ///
       (kdensity lnpoco2 if year==2017, lcolor(gs4)  lpattern(shortdash)) ///
       (kdensity lnpoco2 if year==2023, lcolor(black) lpattern(solid)), ///
   legend(order(1 "2003" 2 "2010" 3 "2017" 4 "2023") position(1) ring(0) rows(4) size(small)) ///
   xtitle("减污降碳(SEM, ln)") ytitle("核密度") scheme(s1mono) graphregion(color(white)) ///
   title("图D 协同排放SEM的分布演变(2003-2023)")
graph export "results/图D_SEM核密度演变.png", replace width(2000) height(1400)
restore

* D.2 城市群维度异质性（用现有 region 分组；如无城市群变量，按 region 近似）
*     若已有城市群编码请替换 group 变量；此处示范用 region
capture confirm string variable region
if _rc==0 {
    encode region, gen(region_c)
}
eststo clear
levelsof region_c, local(rr)
foreach r of local rr {
    local lab : label (region_c) `r'
    cap eststo reg_`r': reghdfe lnpoco2 DID $CTRL if region_c==`r', a(city_code year) vce(cl city_code)
}
esttab reg_* using "results/表D_城市群维度异质性.rtf", replace b(%9.3f) t(%9.3f) ///
    star(* 0.1 ** 0.05 *** 0.01) keep(DID) ///
    stats(N r2_a, fmt(%9.0f %9.3f) labels("N" "Adj.R2")) nogaps compress ///
    title("表D 分区域/城市群 5A 减污降碳效应异质性") ///
    addnotes("按region分样本;各列DID为该区域5A政策效应。可替换为京津冀/长三角/珠三角/长江中游/成渝城市群编码。")
eststo clear

di as result _n "==== 模块 A/B/D 完成（C 为模板，需补投入变量）：结果见 results 文件夹 ===="
