# AI 使用说明

## 工具与运行环境

- AI 工具与模型为 OpenAI Codex 桌面应用中的 GPT-6.2，deepseek v4.1。
- 实际运行环境为 Windows x64、PowerShell、Python 3.13.9（Anaconda）及 CPU。AutoDock Vina 1.2.7 独立程序执行对接，Meeko 0.7.1、RDKit 2025.9.6 和 Gemmi 0.7.4 负责结构处理与结果解析。PyMOL Open-Source 3.2.0a0.post5859 执行光线追踪，Pillow 12.3.0 完成图中文字与图例的排版。

## AI 参与范围

- 对接评估、PyMOL 渲染与图像排版脚本，生成 Vina 配置和项目文档。
- 原始 PDB 与 CCD 文件来自 RCSB。对接姿势和打分来自本地 Vina 运行。RMSD 来自 RDKit 原子映射与 NumPy 坐标计算。三维形状由 PyMOL 根据结构坐标渲染。
- 可视化修订由 AI 调整 PyMOL 视角、结合沟表面、双视角叠合及残基引线。氢键距离由真实坐标测量，图像由deepseek v4.1脚本，导出后渲染与排版。
