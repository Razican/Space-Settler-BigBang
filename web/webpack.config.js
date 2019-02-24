const path = require('path');
const MiniCssExtractPlugin = require("mini-css-extract-plugin");
const envMode = !process.env.NODE_ENV ? 'development' : process.env.NODE_ENV;

module.exports = {
    mode: envMode,
    plugins: [
        new MiniCssExtractPlugin({
            filename: 'css/[name].css',
            chunkFilename: 'css/[id].css',
        })
    ],
    entry: [
        './static/js/script.js',
        './static/css/style.scss',
    ],
    output: {
        publicPath: path.resolve(__dirname, 'dist'),
        path: path.resolve(__dirname, 'dist'),
        filename: 'js/bundle.js'
    },
    module: {
        rules: [{
            test: /\.scss$/,
            use: [{
                loader: MiniCssExtractPlugin.loader,
                options: {
                    sourceMap: envMode == 'development',
                }
            },
            {
                loader: 'css-loader', // translates CSS into CommonJS
                options: {
                    sourceMap: envMode == 'development',
                }
            }, {
                loader: 'sass-loader', // compiles Sass to CSS, using Node Sass by default
                options: {
                    sourceMap: envMode == 'development',
                    outputStyle: envMode == 'development' ? 'expanded' : 'compressed',
                }
            }]
        }],
    },
    devServer: {
        hot: true,
        contentBase: path.join(__dirname, 'dist'),
    },
    devtool: envMode == 'development' ? 'source-map' : '',
};