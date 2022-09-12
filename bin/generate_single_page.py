#!/usr/bin/env python3

import base64
import json
from pathlib import Path


def write_image_data():
    image_data = {}
    images = Path('img/').glob('*.png')
    for image in images:
        with image.open('rb') as f:
            img_data = f.read()
            img_base64 = base64.b64encode(img_data).decode('utf-8')
            image_data[str(image.name)] = f'data:image/png;base64, {img_base64}'

    with Path('js/img_data.js').open('w') as f:
        f.write(f'const img_data = {json.dumps(image_data)}')


def replace_js(file_name: str, dash_html: str) -> str:
    template = f'<script src="{file_name}"></script>'

    with Path(file_name).open() as f:
        js_contents = f.read()

    return dash_html.replace(template, f'<script>\n{js_contents}\n</script>')


def replace_css(file_name: str, dash_html: str) -> str:
    template = f'<link rel="stylesheet" href="{file_name}">'

    with Path(file_name).open() as f:
        css_contents = f.read()

    return dash_html.replace(template, f'<style>\n{css_contents}\n</style>')


def main():
    write_image_data()

    with Path('exp_dash.html').open() as dash_file:
        dash_html = dash_file.read()

    dash_html = replace_js('js/data.js', dash_html)
    dash_html = replace_js('js/log_data.js', dash_html)
    dash_html = replace_js('js/img_data.js', dash_html)
    dash_html = replace_js('js/exp.js', dash_html)
    dash_html = replace_js('js/footer.js', dash_html)

    dash_html = replace_css('style/style.css', dash_html)

    with Path('exp_dash_new.html').open('w') as dash_file:
        dash_file.write(dash_html)


if __name__ == "__main__":
    main()
