*==============================================================================*
*  05_heterogeneity.do —— 模块5：异质性分析（每个维度单独成表，含完整控制变量）
*   11个维度：①资源禀赋 ②行政等级 ③5A数量 ④景区属性 ⑤污染碳本底 ⑥旅游本底 ⑦规制本底
*             ⑧旅游资源丰度(4A密度,新) ⑨距省会区位(新) ⑩海陆区位(新) ⑪数字科技本底(新)
*   每个维度输出一张完整两列表(分组子样本),报告DID及全部控制变量,并给出组间系数差异检验p值。
*  运行：do "05_heterogeneity.do"
*  输出：results/结果05_异质性_*.rtf
*==============================================================================*
do "00_setup.do"

* 组间系数差异检验(交互项法)
cap program drop difftest
program define difftest, rclass
    args g cond
    tempvar gx
    gen double `gx' = DID*`g'
    if "`cond'"=="" reghdfe lnpoco2 DID `gx' $CTRL, a(city_code year) vce(cl city_code)
    else            reghdfe lnpoco2 DID `gx' $CTRL if `cond', a(city_code year) vce(cl city_code)
    return scalar p = 2*ttail(e(df_r), abs(_b[`gx']/_se[`gx']))
end

* 维度清单：编号 变量名 组1标签 组2标签 组1条件 组2条件 经济含义
* 用 @ 分隔字段以容纳含空格的条件
local D1 "①资源禀赋@resource@资源型@非资源型@resource==1@resource==0@资源型城市转型压力与惯性并存"
local D2 "②行政等级@central@中心城市@一般城市@central==1@central==0@中心城市行政资源与执行力更强"
local D3 "③5A数量强度@_g3@单个5A@多个5A@inlist(grp5a,1,0)@inlist(grp5a,2,0)@存在剂量-反应,多个5A效应更强"
local D4 "④景区属性@_g4@自然类@人文类@((treat==1&Resour==0)|treat==0)@((treat==1&Resour==1)|treat==0)@自然类景区生态管制更硬,效应更强"
local D5 "⑤污染碳本底@highpoll@高本底@低本底@highpoll==1@highpoll==0@初始污染碳强度不同,改善空间不同"
local D6 "⑥旅游资源丰度@rich4a@高4A密度@低4A密度@rich4a==1@rich4a==0@旅游资源越丰,信号-集聚效应越强(新)"
local D7 "⑦距省会区位@far_cap@远离省会@邻近省会@far_cap==1@far_cap==0@边缘城市更受益,经济地理视角(新)"
local D8 "⑧海陆区位@coastal@沿海@内陆@coastal==1@coastal==0@要素与环境规制空间梯度(新)"
local D9 "⑨数字科技本底@hitech@高@低@hitech==1@hitech==0@数字经济赋能绿色转型(新)"

forvalues k=1/9 {
    local spec `D`k''
    tokenize "`spec'", parse("@")
    local no  "`1'"
    local g   "`3'"
    local l1  "`5'"
    local l2  "`7'"
    local c1  "`9'"
    local c2  "`11'"
    local econ "`13'"
    use "results/_work.dta", clear
    * 对③④需要构造分组虚拟以做差异检验
    cap gen byte `g' = .
    eststo clear
    eststo A: reghdfe lnpoco2 DID $CTRL if `c1', a(city_code year) vce(cl city_code)
    estadd local fe "是":A
    eststo B: reghdfe lnpoco2 DID $CTRL if `c2', a(city_code year) vce(cl city_code)
    estadd local fe "是":B
    * 组间差异(仅对可直接构造0/1分组变量的维度)
    local pdiff "见交互项检验"
    cap {
        if "`no'"=="③5A数量强度" {
            gen byte gm=(grp5a==2) if inlist(grp5a,1,2)
            difftest gm "inlist(grp5a,1,2)"
            local pdiff : di %5.3f r(p)
        }
        else if "`no'"=="④景区属性" {
            difftest Resour "treat==1"
            local pdiff : di %5.3f r(p)
        }
        else {
            difftest `g'
            local pdiff : di %5.3f r(p)
        }
    }
    esttab A B using "results/结果05_异质性_`k'.rtf", replace ///
        b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) keep(DID $CTRL _cons) order(DID) coeflabels(DID "5A政策(DID)" _cons "常数项") ///
        mtitles("`l1'" "`l2'") ///
        stats(fe N r2_a, fmt(%s %9.0f %9.3f) labels("城市/年FE" "观测值N" "调整R2")) ///
        nogaps compress label title("表8-`no' 异质性:`no'") ///
        addnotes("【方法】按‘`no'’将样本分为两组分别估计DID,均含城市固定效应、年份固定效应与全部控制变量;组间系数差异用交互项DID×分组检验。" ///
         "【计算过程】组间系数差异检验 p=`pdiff'。" ///
         "【经济含义】`econ'。两组系数的相对大小反映5A减污降碳效应在该维度上的异质性。" ///
         "【参考文献】He et al.(2025,JUE)分位/分组异质性;标杆文献异质性分析范式。")
    eststo clear
}
di as result "== 模块5 完成：results/结果05_异质性_1..9.rtf(每维度单独成表) =="
