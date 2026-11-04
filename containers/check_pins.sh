#!/usr/bin/env bash
# Check the published images' pins (containers/pins/*), for the image-pins workflow:
#   1. each recipe input's git tree hash at HEAD equals the one recorded at the
#      commit that built the pinned image (the recipe has not changed since);
#   2. the Snakemake step file references the pinned digest, and no other;
#   3. the pinned manifest is anonymously pullable (HTTP 200).
# Image builds are not byte-reproducible, so this detects "recipe changed since
# the pinned build" rather than rebuilding and comparing. Run from the repo root.
#   bash containers/check_pins.sh [--offline] [pins dir]   (--offline skips 3)
set -uo pipefail
offline=0
if [ "${1:-}" = "--offline" ]; then offline=1; shift; fi
pins=${1:-containers/pins}
fail=0
for f in "$pins"/*; do
    name=$(basename "$f")
    image=$(awk '$1 == "image" {print $2}' "$f")
    digest=$(awk '$1 == "digest" {print $2}' "$f")
    ref=$(awk '$1 == "ref" {print $2}' "$f")
    if [ -z "$image" ] || [ -z "$digest" ] || [ -z "$ref" ]; then
        echo "FAIL $name: image, digest and ref lines are required"; fail=1; continue
    fi
    while read -r _ path want; do
        got=$(git rev-parse "HEAD:$path" 2>/dev/null || echo missing)
        if [ "$got" = "$want" ]; then
            echo "ok   $name: tree $path $got"
        else
            echo "FAIL $name: tree $path is $got, pinned build used $want;" \
                "rebuild the image, then update the digest and $f"
            fail=1
        fi
    done < <(awk '$1 == "tree"' "$f")
    refs=$(grep -o "$image@sha256:[0-9a-f]*" "$ref" | sort -u)
    if [ "$refs" = "$image@$digest" ]; then
        echo "ok   $name: $ref references $digest"
    else
        echo "FAIL $name: $ref references [${refs//$'\n'/ }], expected $image@$digest"
        fail=1
    fi
    if [ "$offline" = 0 ]; then
        repo=${image#ghcr.io/}
        token=$(curl -fsS "https://ghcr.io/token?scope=repository:$repo:pull" |
            python3 -c 'import json, sys; print(json.load(sys.stdin)["token"])' 2>/dev/null)
        if [ -z "$token" ]; then
            echo "FAIL $name: no anonymous pull token for $repo"; fail=1; continue
        fi
        code=$(curl -sS -o /dev/null -w '%{http_code}' -I \
            -H "Authorization: Bearer $token" \
            -H "Accept: application/vnd.oci.image.index.v1+json" \
            -H "Accept: application/vnd.oci.image.manifest.v1+json" \
            -H "Accept: application/vnd.docker.distribution.manifest.list.v2+json" \
            -H "Accept: application/vnd.docker.distribution.manifest.v2+json" \
            "https://ghcr.io/v2/$repo/manifests/$digest")
        if [ "$code" = 200 ]; then
            echo "ok   $name: anonymous manifest HEAD $image@$digest: HTTP $code"
        else
            echo "FAIL $name: anonymous manifest HEAD $image@$digest: HTTP $code"
            fail=1
        fi
    fi
done
[ "$fail" = 0 ] && echo "PINS_OK" || echo "PINS_FAIL"
exit "$fail"
