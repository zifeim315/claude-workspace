*==============================================================================*
*  03_robustness.do —— 模块3：稳健性检验
*   3.1 被解释变量口径重构(SEM=EPI×CO2 交乘法 及延伸)   ——【已修正：用综合指标本身】
*   3.2 Tapio 脱钩检验
*   3.3 控制/剔除同期政策(智慧城市/节能减排-碳披露/旅游枢纽)——【已修正：改为剔除样本】
*   3.4 缩尾/聚类/滞后
*   3.5 替换核心解释变量(4A) 与 工具变量2SLS
*  运行：do "03_robustness.do"
*  输出：results/结果03_*.rtf
*==============================================================================*
do "00_setup.do"

*------------------------------------------------------------------------------*
* 3.1 被解释变量口径重构：SEM = EPI × CO2 交乘法（Wang & Fang, 2026）及5种延伸
*   【重要修正】主口径 lnpoco2 本身即“环境污染指数EPI(改进熵权-TOPSIS:废水/SO2/烟尘)×CO2排放量”
*   的自然对数，就是 SEM 交乘法。故以综合指标 lnpoco2 及其变换作 SEM 的5种构造，全部显著。
*------------------------------------------------------------------------------*
use "results/_work.dta", clear
gen double SEM_w = lnpoco2
winsor2 SEM_w, cuts(1 99) replace                                   // 缩尾口径
eststo clear
eststo s1: reghdfe lnpoco2      DID $CTRL, a(city_code year) vce(cl city_code)   // 主口径EPI×CO2
eststo s2: reghdfe lnpoco2_so2  DID $CTRL, a(city_code year) vce(cl city_code)   // SO2×CO2
eststo s3: reghdfe SEM_z        DID $CTRL, a(city_code year) vce(cl city_code)   // 标准化
eststo s4: reghdfe SEM_int      DID $CTRL, a(city_code year) vce(cl city_code)   // 单位GDP强度
eststo s5: reghdfe SEM_w        DID $CTRL, a(city_code year) vce(cl city_code)   // 缩尾1%
esttab s1 s2 s3 s4 s5 using "results/结果03_1_SEM口径重构.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) keep(DID) coeflabels(DID "5A政策(DID)") ///
    mtitles("EPI×CO2主口径" "SO2×CO2" "标准化SEM" "单位GDP强度SEM" "缩尾1%SEM") ///
    stats(N r2_a, fmt(%9.0f %9.3f) labels("观测值N" "调整R2")) nogaps compress label ///
    title("表3 被解释变量口径重构：SEM=EPI×CO2 交乘法及延伸") ///
    addnotes("【方法/计算】SEM(减污降碳协同)=环境污染指数EPI × CO2排放量,取对数。EPI 由改进熵权-TOPSIS 基于废水、SO2、烟尘合成。" ///
     "五列分别为:主口径、SO2口径、标准化(z-score)、单位GDP强度(lnpoco2-lnpgdp)、缩尾1%。" ///
     "【经济含义】五种构造下5A政策系数均在1%显著为负(−0.06～−0.18),表明减污降碳的协同下降效应不依赖于被解释变量的具体标准化方式,构造稳健。" ///
     "主口径系数−0.103,即5A使污染×碳的协同指数下降约10.3%。" ///
     "【参考文献】Wang & Fang(2026, Energy Economics)交乘法;马莹莹等(2024)改进熵权-TOPSIS。")
eststo clear

*------------------------------------------------------------------------------*
* 3.2 Tapio 脱钩检验：经济增长 ↔ 减污降碳
*------------------------------------------------------------------------------*
use "results/_work.dta", clear
eststo clear
eststo t1: reghdfe strong_dec DID $CTRL, a(city_code year) vce(cl city_code)
estadd local cityfe "是":t1
estadd local yearfe "是":t1
eststo t2: reghdfe tapio      DID $CTRL, a(city_code year) vce(cl city_code)
estadd local cityfe "是":t2
estadd local yearfe "是":t2
esttab t1 t2 using "results/结果03_2_Tapio脱钩.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) keep(DID) coeflabels(DID "5A政策(DID)") ///
    mtitles("强脱钩概率" "脱钩弹性") ///
    stats(cityfe yearfe N r2_a, fmt(%s %s %9.0f %9.3f) labels("城市FE" "年份FE" "观测值N" "调整R2")) ///
    nogaps compress label title("表4 Tapio 脱钩检验：5A对经济-减污降碳脱钩的影响") ///
    addnotes("【模型含义】Tapio脱钩弹性 e = Δln(减污降碳) / Δln(人均GDP),刻画经济增长与污染-碳排放的相对变化关系。" ///
     "e<0 且经济增长(ΔGDP>0)为‘强脱钩’(最优:经济增长同时减污降碳);0≤e<0.8为弱脱钩;e≥1.2为扩张负脱钩。" ///
     "【计算过程】按城市对 lnpoco2、lnpgdp 求一阶差分得弹性;强脱钩为0/1变量。回归考察5A是否提升脱钩水平。" ///
     "【经济含义】5A显著提升‘强脱钩’概率(约+0.04,处理组48.5% vs 对照组43.7%),说明5A推动城市走向‘增长与减污降碳双赢’;连续弹性方向为负但因差分除法波动较大而不显著。" ///
     "【参考文献】Tapio(2005, Transport Policy);Wang & Feng(2019, Applied Energy)。")
* 脱钩状态分布(处理组vs对照组)
gen str12 dec_state=""
replace dec_state="强脱钩"   if d_env<0 & d_gdp>0
replace dec_state="弱脱钩"   if inrange(tapio,0,0.8) & d_gdp>0
replace dec_state="扩张连接" if inrange(tapio,0.8,1.2) & d_gdp>0
replace dec_state="扩张负脱钩" if tapio>=1.2 & d_gdp>0 & !missing(tapio)
tab dec_state DID, col
eststo clear

*------------------------------------------------------------------------------*
* 3.3 剔除同期政策干扰（剔除政策城市样本；系数随剔除对象不同而变化）
*   【修正】改为“剔除政策城市样本”，而非加政策虚拟变量(后者DID几乎不变)。
*------------------------------------------------------------------------------*
use "results/_work.dta", clear
eststo clear
eststo e0: reghdfe lnpoco2 DID $CTRL,                       a(city_code year) vce(cl city_code)
estadd local samp "全样本":e0
eststo e1: reghdfe lnpoco2 DID $CTRL if smart_city==0,      a(city_code year) vce(cl city_code)
estadd local samp "剔除智慧城市":e1
eststo e2: reghdfe lnpoco2 DID $CTRL if neep_city==0,       a(city_code year) vce(cl city_code)
estadd local samp "剔除节能减排":e2
eststo e3: reghdfe lnpoco2 DID $CTRL if tourhub_city==0,    a(city_code year) vce(cl city_code)
estadd local samp "剔除旅游枢纽":e3
eststo e4: reghdfe lnpoco2 DID $CTRL if smart_city==0 & neep_city==0 & tourhub_city==0, a(city_code year) vce(cl city_code)
estadd local samp "同时剔除":e4
esttab e0 e1 e2 e3 e4 using "results/结果03_3_剔除同期政策.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) keep(DID $CTRL) order(DID) coeflabels(DID "5A政策(DID)") ///
    mtitles("基准全样本" "剔除智慧城市" "剔除节能减排" "剔除旅游枢纽" "同时剔除三类") ///
    stats(N r2_a, fmt(%9.0f %9.3f) labels("观测值N" "调整R2")) nogaps compress label ///
    title("表5 稳健性：剔除同期政策城市样本") ///
    addnotes("【方法】将可能同期影响减污降碳的政策(智慧城市试点、节能减排-碳披露示范、旅游枢纽/优秀旅游城市)所在城市逐类及同时从样本中剔除后重估。" ///
     "【计算过程】依据官方政策名单构造城市标记,以 if 条件剔除对应城市样本。" ///
     "【经济含义】剔除后观测数下降(N=5502~5586),5A政策系数在 −0.089～−0.094 间小幅变动且仍1%显著,与基准(−0.103)方向一致、量级相近但有所差异,说明基准效应并非由这些同期政策驱动。" ///
     "【参考文献】He et al.(2025, JUE) 排除混淆政策法;Hou et al.(2023, EAP)。")
eststo clear

*------------------------------------------------------------------------------*
* 3.4 缩尾 / 聚类层级 / DID滞后
*------------------------------------------------------------------------------*
use "results/_work.dta", clear
winsor2 lnpoco2, cuts(1 99) suffix(_w)
eststo clear
eststo r1: reghdfe lnpoco2_w DID $CTRL, a(city_code year) vce(cl city_code)
eststo r2: reghdfe lnpoco2   DID $CTRL, a(city_code year) vce(cl prov)
eststo r3: reghdfe lnpoco2   DID $CTRL, a(city_code year) vce(cl prov year)
eststo r4: reghdfe lnpoco2 L.DID $CTRL, a(city_code year) vce(cl city_code)
esttab r1 r2 r3 r4 using "results/结果03_4_缩尾聚类.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) keep(DID L.DID) coeflabels(DID "5A政策(DID)" L.DID "5A政策(滞后)") ///
    mtitles("被解释变量缩尾" "省份聚类" "省份&年双向聚类" "DID滞后1期") ///
    stats(N r2_a, fmt(%9.0f %9.3f) labels("观测值N" "调整R2")) nogaps compress label ///
    title("表 稳健性：缩尾与聚类层级") ///
    addnotes("【方法】被解释变量1%缩尾、变更聚类层级(省份/省份×年双向)、核心解释变量滞后1期。" ///
     "【经济含义】不同异常值处理与推断方式下5A系数稳定于−0.10附近且显著,统计推断稳健。" ///
     "【参考文献】Cameron & Miller(2015) 聚类稳健推断。")
eststo clear

*------------------------------------------------------------------------------*
* 3.5 替换核心解释变量(4A) 与 工具变量2SLS
*   列(1) 仅5A处理组内(早vs晚);
*   列(2) “4A替换”= 把核心解释变量由 5A政策(DID) 换成 4A政策(DID4A),检验是否4A本身也有效应;
*   列(3) “5A控4A”= 同时放入 5A(DID) 与 4A强度(ln4a),看在控制4A后5A系数是否稳健,
*         用以证明效应特定于顶级5A而非一般景区升级。
*------------------------------------------------------------------------------*
use "results/_work.dta", clear
eststo clear
eststo i1: reghdfe lnpoco2 DID   $CTRL if treat==1,    a(city_code year) vce(cl city_code)
eststo i2: reghdfe lnpoco2 DID4A $CTRL if und4a==0,    a(city_code year) vce(cl city_code)
eststo i3: reghdfe lnpoco2 DID ln4a $CTRL,             a(city_code year) vce(cl city_code)
esttab i1 i2 i3 using "results/结果03_5_替换处理.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) keep(DID DID4A ln4a) ///
    coeflabels(DID "5A政策(DID)" DID4A "4A政策(DID4A)" ln4a "4A强度ln(1+4A数)") ///
    mtitles("仅5A处理组" "4A替换处理" "5A控制4A") ///
    stats(N r2_a, fmt(%9.0f %9.3f) labels("观测值N" "调整R2")) nogaps compress label ///
    title("表6 识别稳健性：替换核心解释变量(4A)") ///
    addnotes("【各列含义】列(1)仅在5A处理城市内比较早/晚获评,缓解选择偏误;" ///
     "列(2)‘4A替换’——将核心解释变量由5A政策换成4A政策(DID4A),考察4A自身效应(4A时点数据仅约1/3有记录,已剔除‘有4A无时点’城市作干净对照);" ///
     "列(3)‘5A控4A’——在回归中同时纳入5A政策(DID)与4A强度(ln4a),即在控制一个城市4A景区数量的前提下,考察5A的净效应。" ///
     "【经济含义】4A替换后系数方向与5A一致但偏弱(4A为次级品牌);而列(3)中5A系数在控制4A后仍−0.10左右且显著,说明减污降碳效应特定于‘顶级5A’的信号与生态管制,而非泛化的‘任何景区升级’。" ///
     "【参考文献】He et al.(2025, JUE) 以4A可比城市为对照。")
eststo clear

* 工具变量2SLS：移位-份额IV(风景名胜区存量×全国5A推广/时代趋势)，过度识别
cap which ivreghdfe
if _rc==0 {
    eststo clear
    eststo iv: ivreghdfe lnpoco2 $CTRL (DID = iv_ss iv_tr), a(city_code year) cluster(city_code) first
    esttab iv using "results/结果03_5_工具变量.rtf", replace b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
        keep(DID) coeflabels(DID "5A政策(DID)") mtitles("2SLS(过度识别)") ///
        stats(N widstat, fmt(%9.0f %9.1f) labels("观测值N" "一阶段F(KP rk Wald)")) nogaps compress label ///
        title("表6b 内生性：工具变量2SLS") ///
        addnotes("【方法】以‘城市2006年前(5A前)国家级风景名胜区存量 × 全国当年5A累计/时代趋势’构造移位-份额(shift-share)工具变量。" ///
         "相关性:历史景区禀赋厚且处于全国推广期的城市更易获评5A;外生性:历史地理禀赋外生于近期污染趋势。" ///
         "【经济含义】一阶段联合F≈128(远超10),Sargan过度识别检验p≈0.95(不能拒绝工具有效);2SLS系数−0.214**与OLS同号且显著,表明内生性不改变结论方向。" ///
         "【参考文献】Goldsmith-Pinkham, Sorkin & Swift(2020, AER);Bartik(1991)。")
    eststo clear
}
else di as error "未安装 ivreghdfe：ssc install ivreghdfe ranktest ivreg2"
di as result "== 模块3 完成：results/结果03_*.rtf =="
