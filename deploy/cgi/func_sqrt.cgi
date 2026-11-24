#!/home/cf326296/miniconda3/bin/python3
"""func_sqrt API の CGI エントリポイント。

GET  /func_sqrt?f=sinh(x)&n=11
POST /func_sqrt  (application/x-www-form-urlencoded または application/json)
"""
import json
import os
import signal
import sys
from urllib.parse import parse_qs

# パッケージの場所 (公開ディレクトリの外に置く)
sys.path.insert(0, os.path.expanduser("~/apps/func_sqrt"))

TIME_LIMIT = 10  # 秒


def respond(status, body, cache=False):
    reasons = {200: "OK", 400: "Bad Request", 405: "Method Not Allowed",
               500: "Internal Server Error", 504: "Gateway Timeout"}
    out = sys.stdout
    out.write(f"Status: {status} {reasons.get(status, '')}\r\n")
    out.write("Content-Type: application/json; charset=utf-8\r\n")
    out.write("Access-Control-Allow-Origin: *\r\n")
    out.write("Access-Control-Allow-Methods: GET, POST, OPTIONS\r\n")
    out.write("Access-Control-Allow-Headers: Content-Type\r\n")
    out.write("Cache-Control: " + ("public, max-age=86400" if cache else "no-store") + "\r\n")
    out.write("\r\n")
    out.flush()
    sys.stdout.buffer.write(json.dumps(body, ensure_ascii=False, indent=2).encode("utf-8"))


def read_params():
    params = {k: v[0] for k, v in parse_qs(os.environ.get("QUERY_STRING", "")).items()}
    if os.environ.get("REQUEST_METHOD") == "POST":
        length = min(int(os.environ.get("CONTENT_LENGTH") or 0), 10000)
        raw = sys.stdin.buffer.read(length).decode("utf-8", "replace")
        if "json" in os.environ.get("CONTENT_TYPE", ""):
            data = json.loads(raw or "{}")
            if not isinstance(data, dict):
                raise ValueError("JSON はオブジェクトで送ってください")
            params.update({k: str(v) for k, v in data.items()})
        else:
            params.update({k: v[0] for k, v in parse_qs(raw).items()})
    return params


def on_timeout(signum, frame):
    raise TimeoutError


def main():
    method = os.environ.get("REQUEST_METHOD", "GET")
    if method == "OPTIONS":
        return respond(200, {"ok": True})
    if method not in ("GET", "POST"):
        return respond(405, {"ok": False, "error": "GET か POST で呼び出してください"})

    from func_sqrt.web import handle, usage

    try:
        params = read_params()
    except ValueError as exc:
        return respond(400, {"ok": False, "error": str(exc)})
    if not any(params.get(k) for k in ("f", "derivs", "coeffs")):
        return respond(200, usage(), cache=True)

    signal.signal(signal.SIGALRM, on_timeout)
    signal.alarm(TIME_LIMIT)
    try:
        status, body = handle(params)
    except TimeoutError:
        status, body = 504, {"ok": False, "error": f"計算が {TIME_LIMIT} 秒以内に終わりませんでした。次数を下げてください。"}
    except Exception:
        status, body = 500, {"ok": False, "error": "サーバー内部でエラーが発生しました"}
    finally:
        signal.alarm(0)
    respond(status, body, cache=(status == 200))


if __name__ == "__main__":
    main()
