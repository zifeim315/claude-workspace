*==============================================================================*
*  00_setup.do —— 公共数据准备（所有分模块 do 文件在开头 do 本文件即可）
*  数据：data_spatial_v2.dta（289城市×2003-2023=6069）
*  作用：载入数据、构造所有派生变量与全局宏、保存 results/_work.dta 供各模块复用
*  用法：单独运行本文件，或在各分模块顶部  do "00_setup.do"
*  外部命令：ssc install reghdfe ftools estout winsor2 boottest ivreghdfe ranktest ivreg2 ///
*            csdid drdid bacondecomp psmatch2 xsmle
*==============================================================================*
clear all
set more off
set matsize 2000
cap mkdir results
global MASTERDIR "results"

use "data_spatial_v2.dta", clear
xtset city_code year

*--- 控制变量与固定效应识别单元 ---*
cap gen double lnfin = ln(fin)
label var lnfin "金融发展水平(存贷余额/GDP,对数)"
global CTRL   "lnpgdp lndensity urban struc2 gov tech human lnfin"   // 全控制(8个)
global CTRLns "lnpgdp lndensity urban gov tech human lnfin"          // 剔除结构变量(用于结构类机制)
cap egen prov = group(province)
cap egen reg  = group(region)
cap egen provyear = group(prov year)
cap egen regyear  = group(reg  year)

*--- 被解释变量口径（SEM=EPI×CO2 的多种构造；主口径 lnpoco2 即熵权-TOPSIS 交乘）---*
* 说明：lnpoco2 = ln(环境污染指数EPI[改进熵权-TOPSIS:废水/SO2/烟尘] × CO2排放量)，即 SEM 交乘法主口径
qui sum lnpoco2
gen double SEM_z   = (lnpoco2 - r(mean))/r(sd)                       // 标准化口径
gen double SEM_int = lnpoco2 - lnpgdp                                // 单位GDP强度口径
label var SEM_z   "减污降碳(标准化SEM)"
label var SEM_int "减污降碳(单位GDP强度SEM)"

*--- 公众环境关注度(中介5)：百度指数标准化 ---*
qui sum haze
cap gen double zhaze = (haze - r(mean))/r(sd)
qui sum search
cap gen double zatt  = (search - r(mean))/r(sd)
label var zhaze "公众雾霾关注度(标准化)"
label var zatt  "公众环境关注度(标准化)"

*--- 工具变量：风景名胜区存量 × 全国5A推广/时代趋势（移位-份额）---*
cap bysort year: egen nat5a = total(DID)
gen double iv_ss = scenic_pre * nat5a
gen double iv_tr = scenic_pre * max(year-2006,0)
label var iv_ss "IV1:风景名胜区存量×全国5A推广"
label var iv_tr "IV2:风景名胜区存量×(year-2006)"

*--- 同期政策“剔除样本”标记（智慧城市/节能减排-碳披露/旅游枢纽）---*
* 名单以官方文件为准，可在此增删；生效年近似取政策首批年。
gen byte smart_city = 0
foreach c in 北京市 上海市 南京市 杭州市 宁波市 温州市 无锡市 扬州市 武汉市 广州市 深圳市 成都市 西安市 ///
             厦门市 青岛市 济南市 郑州市 合肥市 福州市 昆明市 太原市 呼和浩特市 沈阳市 大连市 哈尔滨市 秦皇岛市 石家庄市 {
    replace smart_city=1 if city=="`c'"
}
gen byte neep_city = 0
foreach c in 北京市 深圳市 重庆市 杭州市 长沙市 贵阳市 吉林市 新余市 石家庄市 唐山市 铁岭市 齐齐哈尔市 ///
             南京市 南昌市 临沂市 鹤壁市 柳州市 广元市 六盘水市 韶关市 宁波市 铜陵市 厦门市 {
    replace neep_city=1 if city=="`c'"
}
gen byte tourhub_city = 0
foreach c in 北京市 天津市 上海市 重庆市 杭州市 苏州市 南京市 成都市 西安市 桂林市 三亚市 厦门市 青岛市 ///
             大连市 昆明市 丽江市 张家界市 黄山市 承德市 秦皇岛市 洛阳市 开封市 九江市 泰安市 宜昌市 {
    replace tourhub_city=1 if city=="`c'"
}
label var smart_city   "智慧城市试点城市"
label var neep_city    "节能减排(碳披露)示范城市"
label var tourhub_city "旅游枢纽/优秀旅游城市"

*--- Tapio 脱钩变量 ---*
gen double d_env = D.lnpoco2
gen double d_gdp = D.lnpgdp
gen double tapio = d_env/d_gdp
winsor2 tapio, cuts(2.5 97.5) replace
gen byte strong_dec = (d_env<0 & d_gdp>0) if !missing(d_env,d_gdp)
label var strong_dec "强脱钩(经济增长且减污降碳)"
label var tapio      "Tapio脱钩弹性"

*--- 异质性新维度分组基元 ---*
gen byte coastal = strpos(region,"东")>0
bysort city_code: egen _dc = mean(dist_cap)
qui sum _dc, detail
gen byte far_cap = _dc > r(p50)
bysort city_code: egen _t4 = max(num4a_cum)
qui sum _t4, detail
gen byte rich4a = _t4 > r(p50)
bysort city_code: egen _pt = mean(cond(DID==0, tech, .))
qui sum _pt, detail
gen byte hitech = _pt > r(p50)
label var coastal "沿海(东部)"
label var far_cap "距省会较远"
label var rich4a  "旅游资源丰度(4A密度高)"
label var hitech  "数字/科技本底高"

save "results/_work.dta", replace
di as result "== 00_setup 完成：results/_work.dta 已就绪，全局宏 CTRL/CTRLns 已设定 =="
