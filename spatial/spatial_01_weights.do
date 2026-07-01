*==============================================================================*
*  空间溢出分析 · 步骤1：读入空间权重矩阵
*  说明：所有权重矩阵已在 Python 端按 289 个城市、按 city_code 升序构造并行标准化，
*        导出为 CSV（无表头，289×289）。本 do 将其读入 Stata 矩阵，供 xsmle 使用。
*  运行前需安装：ssc install xsmle ;  ssc install xthreg ;  ssc install spmat
*  目录约定：本文件夹内含 data_spatial.dta 、各 W_*.csv 、city_order.csv
*==============================================================================*
clear all
set more off
* 如有需要，把工作目录切到本文件夹：
* cd "你的路径/spatial"

*--- 定义一个把 CSV 读成 Stata 矩阵的小程序 ---*
capture program drop loadW
program define loadW
    args csv matname
    preserve
        import delimited using "`csv'", clear varnames(nonames)
        mkmat _all, matrix(`matname')
    restore
end

*--- 逐个读入（矩阵已行标准化，无需再 normalize）---*
loadW "W_adj.csv"      Wadj        // ① 邻接矩阵
loadW "W_geo.csv"      Wgeo        // ② 地理距离(反距离平方)
loadW "W_econ.csv"     Wecon       // ③ 经济距离(人均GDP差)
loadW "W_econgeo.csv"  Wecongeo    // ④ 经济地理嵌套矩阵(主推)
* 距离分段矩阵（溢出范围检验用）
foreach km in 150 200 250 300 350 400 450 500 {
    loadW "W_d`km'.csv" Wd`km'
}
di as result "== 全部空间权重矩阵已读入内存（Wadj Wgeo Wecon Wecongeo Wd150...Wd500）=="
