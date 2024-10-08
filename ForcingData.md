# 驱动数据

## 气象驱动

- `Rn`: Net radiation, [W/m-2]
- `Ta`: air temperature, [degC]
- `Pe`: precipitation, [mm/d]
- `Pa`: air pressure, [W/m-2]

## 计算得到的变量

- `Tas`: 有效积温，土壤模块

- `G`  : 土壤热通量, [W m-2]

    ```r
    Gi = 0.4 * Rn * exp(-0.5 .* LAI); # Choudhury et al., 1987
    ```

## 植被驱动

- `LAI`: leaf area index
- `VOD`: vegetation optical depth, (0-1)

## 中间变量

- `wa` : 三层土壤含水量
- `zwg`: groundwater table depth, mm
- `snp`: snow package, mm


# 可以替换的变量

## 土壤参数

除了`theta_c`: http://globalchange.bnu.edu.cn/research/soil2

## 根系参数

Hydrologic Regulation of Plant Rooting Depth
