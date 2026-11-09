// FastNLR desktop design tokens.
// Author: Jiwen Zhao (https://github.com/CropCoder)
/** @type {import('tailwindcss').Config} */
export default {
  content: ["./index.html", "./src/**/*.{ts,tsx}"],
  theme: {
    extend: {
      colors: {
        mint: "#5cc5b2",
        sun: "#eccc68",
        iris: "#9c88ff",
        coral: "#ff6b6b",
        ink: "#112031",
        canvas: "#f4f7f6"
      },
      fontFamily: {
        sans: ["Inter", "Segoe UI", "Arial", "sans-serif"],
        chart: ["Arial", "sans-serif"]
      }
    }
  },
  plugins: []
};
