## overview
仅有一个类别型变量,称为单因素方差分析(one-way ANOVA),或单因素组间方差分析。方差分析主要通过F检验来进行效果评测,若F检验显著,
则说明group均值不同。

## function
| design        | expression   |
| --------    | :----:  |
| 单因素ANOVA                                |        y ~ A |
| 含单个协变量的单因素ANCOVA                    |    y ~ x + A|
| 双因素ANOVA                                |       y ~ A * B |
| 含两个协变量的双因素ANCOVA                   |     y ~ x1 + x2 + A*B |
| 随机化区组                                   |       y ~ B + A ( B 是区组因子) |
| 单因素组内ANOVA                              |       y ~ A + Error(Subject/A) |
| 含单个组内因子(W)和单个组间因子(B)的重复测量ANOVA |    y ~ B * W + Error(Subject/W) |
