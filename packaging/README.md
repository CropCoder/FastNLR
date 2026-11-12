# FastNLR TBtools 插件打包

`Build-TBtools-Plugin-V1.sh` 把 `fastnlr-server` 打包为 TBtools 可拉起的本地 Web 插件：
打开时在本机启动服务并跳转到 Web UI，离线可用。

用法：

```bash
bash packaging/Build-TBtools-Plugin-V1.sh macos v1.9.2
bash packaging/Build-TBtools-Plugin-V1.sh windows v1.9.2
```

产物为 `dist/fastnlr-server-tbtools-<target>.plugin.zip`。正式提交 TBtools 插件商店的
`.plugin` 包请用 TBtools 的 Plugin Creator 包装本目录的二进制与 launch 脚本。
