import re
import json
import sys

# Download the html tag list from https://www.dicomlibrary.com/dicom/dicom-tags/
default_html_path = '/Users/kliron/Downloads/DICOMlib.html'
default_out = '../../dicom_tags.json'


def get_tags(path: str, out: str) -> None:
    dicom_tags = {}
    with open(path, 'r') as f:
        lines = f.readlines()
        rgx = re.compile(r'\s*<tr><td.*?>\(([0-9A-Za-z]{4}),([0-9A-Za-z]{4})\)</td><td.*?>.*?</td><td.*?>(.*?)<')
        for line in lines:
            m = rgx.match(line)
            if not m:
                if line and not line.isspace():
                    print(f'Stopping after no match at line:\n {line}')
                    break
                continue
            dicom_tags[m.group(1) + '|' + m.group(2)] = m.group(3)

    with open(out, 'w') as w:
        w.write(json.dumps(dicom_tags, sort_keys=True, indent=4))


if __name__ == '__main__':
    if len(sys.argv) == 3:
        r_file = sys.argv[1]
        w_file = sys.argv[2]
    elif len(sys.argv) == 2:
        r_file = sys.argv[1]
        w_file = default_out
    else:
        r_file = default_html_path
        w_file = default_out

    get_tags(r_file, w_file)

