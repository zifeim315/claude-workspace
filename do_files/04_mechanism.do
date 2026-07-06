*==============================================================================*
*  04_mechanism.do —— 模块4：作用机制（三步中介法 + Sobel）
*   五条机制：①产业结构高级化(三产集聚) ②能源强度优化 ③绿色技术创新 ④环境规制强化 ⑤公众环境关注度
*   每条机制输出“三步法”三列表(对标0630文档表范式)：第一步 lnpoco2、第二步 中介、第三步 lnpoco2。
*  运行：do "04_mechanism.do"
*  输出：results/结果04_机制_*.rtf
*==============================================================================*
do "00_setup.do"

* 机制清单：标签 中介变量 控制集(CF=全控制;CA=剔除结构变量以避免坏控制,江艇2022)
* 结构类机制(三产集聚)用 CA；其余用 CF。
local MECH `" "三产集聚 ter_gdp CA H2a" "能源强度 lnelec_gdp CF H2b" "绿色创新 lnpatapp CF H2c" "环境规制 er CF H2d" "公众关注 zhaze CF H2e" "'

di as result _n "{hline 92}"
di as text %-12s "机制" %10s "a(DID→M)" %12s "b(M→Y)" %11s "c'(DID)" %11s "a*b" %10s "Sobel z" %8s "p"
di as text "{hline 92}"

foreach row of local MECH {
    gettoken tag row : row
    gettoken med row : row
    gettoken cid row : row
    gettoken hyp row : row
    if "`cid'"=="CF" local C "$CTRL"
    else             local C "$CTRLns"
    use "results/_work.dta", clear
    eststo clear
    eststo m1: reghdfe lnpoco2 DID `C',       a(city_code year) vce(cl city_code)
    scalar cco=_b[DID]
    eststo m2: reghdfe `med'   DID `C',       a(city_code year) vce(cl city_code)
    scalar aco=_b[DID]
    scalar sea=_se[DID]
    eststo m3: reghdfe lnpoco2 DID `med' `C', a(city_code year) vce(cl city_code)
    scalar bco=_b[`med']
    scalar seb=_se[`med']
    scalar cpr=_b[DID]
    scalar indd=aco*bco
    scalar zso=indd/sqrt(bco^2*sea^2+aco^2*seb^2)
    scalar pso=2*(1-normal(abs(zso)))
    local sg=cond(abs(zso)>2.58,"***",cond(abs(zso)>1.96,"**",cond(abs(zso)>1.65,"*","ns")))
    esttab m1 m2 m3 using "results/结果04_机制_`tag'.rtf", replace ///
        b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) keep(DID `med' _cons) ///
        coeflabels(DID "5A政策(DID)" _cons "常数项") ///
        mtitles("第一步 lnpoco2" "第二步 `tag'" "第三步 lnpoco2") ///
        stats(N r2_a, fmt(%9.0f %9.3f) labels("观测值N" "调整R2")) nogaps compress label ///
        title("表7 机制检验(`hyp')：`tag' 三步法") ///
        addnotes("【方法】温忠麟式三步中介法:第一步 Y=cDID(总效应);第二步 M=aDID;第三步 Y=c'DID+bM。间接效应=a×b,以Sobel检验其显著性。" ///
         "结构类机制控制集剔除struc2以避免‘坏控制’(江艇,2022)。" ///
         "【计算过程】a=`=string(aco,"%6.3f")',b=`=string(bco,"%6.3f")',间接效应a*b=`=string(indd,"%6.4f")',Sobel z=`=string(zso,"%5.2f")',p=`=string(pso,"%5.3f")'(`sg')。" ///
         "【经济含义】5A通过‘`tag'’渠道影响减污降碳:该中介对DID显著响应(a路径),且中介对减污降碳有显著作用(b路径),间接效应显著,机制`hyp'成立。" ///
         "【参考文献】温忠麟等(2004);江艇(2022,中国工业经济);He et al.(2025,JUE)部门溢出。")
    di as text %-12s "`tag'" as result %10.3f aco %12.3f bco %11.3f cpr %11.4f indd %10.3f zso %8.3f pso as text " `sg'"
    eststo clear
}
di as text "{hline 92}"

*------------------------------------------------------------------------------*
* 4.2 Bootstrap 间接效应（百分位95%置信区间，500次，城市聚类重抽样）
*   稳健于中介效应非正态；95%CI不含0 即间接效应显著。
*------------------------------------------------------------------------------*
*   注：城市聚类重抽样会使同一城市被抽多次而产生“重复面板ID”，故用 idcluster(newid)
*       为每次抽到的城市赋唯一新ID,并在 reghdfe 中吸收 newid(而非 city_code),避免 r(451)。
cap program drop bootmed
program bootmed, rclass
    reghdfe ${med} DID ${cc}, a(newid year)
    local a = _b[DID]
    reghdfe lnpoco2 DID ${med} ${cc}, a(newid year)
    return scalar ind = `a'*_b[${med}]
end
tempname BM
postfile `BM' str16 mech double(ind lo hi) using "results/_bootmed.dta", replace
di as result _n "{hline 74}"
di as text %-14s "机制(中介)" %12s "间接效应a*b" %14s "Boot 95%CI下" %14s "Boot 95%CI上" "  显著"
di as text "{hline 74}"
local bi = 0
foreach row of local MECH {
    gettoken tag row : row
    gettoken med row : row
    gettoken cid row : row
    gettoken hyp row : row
    local ++bi
    if "`cid'"=="CF" global cc "$CTRL"
    else            global cc "$CTRLns"
    global med "`med'"
    use "results/_work.dta", clear
    qui bootstrap ind=r(ind), reps(500) seed(20250624) cluster(city_code) idcluster(newid) ///
        saving("results/_br`bi'.dta", replace) nodots: bootmed
    scalar pe = _b[ind]
    preserve
        use "results/_br`bi'.dta", clear          // 自抽样重复值,手工求百分位CI(稳健)
        _pctile ind, p(2.5 97.5)
        scalar lo = r(r1)
        scalar hi = r(r2)
    restore
    local sig = cond(lo*hi>0,"显著","不显著")
    post `BM' ("`tag'") (pe) (lo) (hi)
    di as text %-14s "`tag'" as result %12.4f pe %14.4f lo %14.4f hi as text "  `sig'"
}
postclose `BM'
di as text "{hline 74}"
* 导出 Bootstrap 结果表
preserve
    use "results/_bootmed.dta", clear
    format ind lo hi %9.4f
    gen CI = "[" + string(lo,"%6.4f") + ", " + string(hi,"%6.4f") + "]"
    gen sig = cond(lo*hi>0,"显著(不含0)","不显著")
    list mech ind CI sig, sep(0) noobs
    export excel mech ind lo hi CI sig using "results/结果04_机制Bootstrap.xlsx", replace first(var)
restore

di as result "== 模块4 完成：results/结果04_机制_*.rtf 与 结果04_机制Bootstrap.xlsx =="
