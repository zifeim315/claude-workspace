*==============================================================================*
*  run_all.do —— 主控文件：一次性运行全部模块（也可分别单独运行各模块以节省时间）
*  用法：先 cd 到本文件夹(内含 data_spatial_v2.dta 与各模块 do)，再  do "run_all.do"
*  提示：各模块均可独立运行(每个模块开头会自动 do "00_setup.do")。
*        空间模块(06)需安装 xsmle，耗时较长，可单独运行。
*==============================================================================*
do "00_setup.do"          // 公共数据准备 -> results/_work.dta
do "01_descriptive.do"    // 表1 描述性统计
do "02_baseline.do"       // 表2 基准FE阶梯 + 平行趋势 + 安慰剂 + PSM-DID
do "03_robustness.do"     // 表3 SEM口径 / 表4 Tapio脱钩 / 表5 剔除政策 / 表6 4A与IV
do "04_mechanism.do"      // 表7 五条机制三步法
do "05_heterogeneity.do"  // 表8 十一维异质性(每维一表)
do "06_spatial.do"        // 表9 LISA + SDM/Moran/距离衰减/空间DID
di as result _n "================ 全部模块运行完成，结果见 results/ 文件夹 ================"
