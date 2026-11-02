# FastNLR 桌面版（Tauri v2 脚手架）

本目录是 FastNLR 桌面版的 Tauri v2 脚手架，复用 `nlr-service` 作为执行层。
当前为代码脚手架，未在 macOS/Windows 上完成打包验证。

## 构建

```bash
cd app
npm install
npm run tauri build
```

前置条件：Node.js 20+、Rust 1.77+、以及对应平台的系统依赖
（macOS 无额外依赖，Windows 需要 WebView2，Linux 需要 webkit2gtk-4.1）。

## 结构

- `index.html`：极简前端占位，正式版将替换为与 Web 服务共用的 React 界面。
- `src-tauri/`：Rust 后端，通过 `nlr-service` 执行基因组/蛋白注释。
