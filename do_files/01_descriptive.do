*==============================================================================*
*  01_descriptive.do —— 模块1：描述性统计
*  运行：do "01_descriptive.do"（内部自动 do 00_setup.do）
*  输出：results/结果01_描述性统计.rtf
*==============================================================================*
do "00_setup.do"
use "results/_work.dta", clear

eststo clear
estpost summarize lnpoco2 DID $CTRL er ter_gdp sec_gdp lnelec_gdp lnpatapp zhaze, detail
esttab using "results/结果01_描述性统计.rtf", replace ///
    cells("count(fmt(0)) mean(fmt(3)) sd(fmt(3)) min(fmt(3)) p50(fmt(3)) max(fmt(3))") ///
    noobs nonumber label collabels("N" "均值" "标准差" "最小值" "中位数" "最大值") ///
    title("表1 主要变量描述性统计") ///
    addnotes("【方法】全样本 289 个地级市 × 2003–2023 年，N=6069。" ///
             "【计算过程】被解释变量 lnpoco2 = ln(环境污染指数 EPI × CO2 排放量)，其中 EPI 由改进熵权–TOPSIS 法基于废水、SO2、烟尘三项合成。" ///
             "【经济含义】lnpoco2 均值 4.113、标准差 1.797，城市间减污降碳水平差异较大；DID 均值 0.331，样本期约1/3的城市-年份处于5A政策实施状态。" ///
             "【参考文献】马莹莹等(2024);马彦瑞等(2024)。")
eststo clear
di as result "== 模块1 完成：results/结果01_描述性统计.rtf =="
