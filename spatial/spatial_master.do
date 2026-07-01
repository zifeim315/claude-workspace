*==============================================================================*
*  空间溢出分析 · 主控文件（一键复现）
*  作者：（研究者）  数据：289 个地级市 × 2003–2023（完全平衡面板）
*  被解释变量 lnpoco2（减污降碳，越小越好）；核心 DID（5A 景区政策）。
*
*  运行前准备：
*   1) 安装用户命令： ssc install xsmle ;  ssc install xthreg ;  ssc install spmat
*   2) 把工作目录切到本文件夹（内含 data_spatial.dta 与各 W_*.csv）：
*        cd "你的路径/spatial"
*   3) 依次运行本主控即可。
*
*  权重矩阵（289×289，已按 city_code 升序、行标准化）：
*   Wadj 邻接 | Wgeo 地理距离 | Wecon 经济距离 | Wecongeo 经济地理嵌套(主推)
*   Wd150…Wd500 距离带（溢出范围）
*==============================================================================*
clear all
set more off

do "spatial_01_weights.do"     // 读入全部空间权重矩阵
do "spatial_02_sdm.do"         // SDM 选矩阵 + 直接/间接/总效应
do "spatial_03_distance.do"    // 距离衰减：溢出范围
do "spatial_04_threshold.do"   // 门槛/分段异质性溢出
do "spatial_05_mechanism.do"   // 溢出传导机制（产业结构邻近性）

di as result _n "================ 空间溢出全部分析完成，结果见 results_*.rtf / *.png ================"
