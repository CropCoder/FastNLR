// FastNLR desktop test-runner configuration.
// Author: Jiwen Zhao (https://github.com/CropCoder)
import { defineConfig } from "vitest/config";

export default defineConfig({
  test: {
    environment: "jsdom",
    setupFiles: "./src/test/setup.ts",
    css: true,
  },
});
