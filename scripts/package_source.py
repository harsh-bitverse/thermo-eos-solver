import os
import sys
import zipfile
import subprocess
import shutil

EXCLUDE_DIRS = {
    'build', '.git', 'dist', '__pycache__', '.vscode', '.idea', 'tmp_test_extract'
}
EXCLUDE_EXTS = {
    '.exe', '.obj', '.o', '.a', '.dll', '.pdb', '.pyc', '.ninja', '.log'
}

def create_source_archive(output_zip="dist/thermo-eos-solver-source.zip"):
    root_dir = os.path.abspath(".")
    os.makedirs(os.path.dirname(output_zip), exist_ok=True)

    print(f"[*] Packaging clean source archive to {output_zip}...")
    files_added = 0
    with zipfile.ZipFile(output_zip, 'w', zipfile.ZIP_DEFLATED) as zf:
        for dirpath, dirnames, filenames in os.walk(root_dir):
            # Prune excluded directories
            dirnames[:] = [d for d in dirnames if d not in EXCLUDE_DIRS]

            for filename in filenames:
                _, ext = os.path.splitext(filename)
                if ext.lower() in EXCLUDE_EXTS:
                    continue
                if filename.endswith(".zip"):
                    continue

                full_path = os.path.join(dirpath, filename)
                rel_path = os.path.relpath(full_path, root_dir)

                # Store inside a top-level folder in the zip
                archive_name = os.path.join("thermo-eos-solver", rel_path).replace("\\", "/")
                zf.write(full_path, archive_name)
                files_added += 1

    print(f"[SUCCESS] Packaged {files_added} files into {output_zip} (size: {os.path.getsize(output_zip):,} bytes)")
    return output_zip

def verify_extracted_build(zip_path):
    temp_dir = os.path.abspath("tmp_test_extract")
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir, ignore_errors=True)
    os.makedirs(temp_dir, exist_ok=True)

    print(f"[*] Verifying archive extraction into temporary directory: {temp_dir}...")
    with zipfile.ZipFile(zip_path, 'r') as zf:
        zf.extractall(temp_dir)

    extracted_root = os.path.join(temp_dir, "thermo-eos-solver")
    print(f"[*] Configuring CMake in extracted directory: {extracted_root}...")

    # Configure
    cmd_cfg = ["cmake", "-B", "build", "-S", ".", "-DCMAKE_BUILD_TYPE=Release"]
    res_cfg = subprocess.run(cmd_cfg, cwd=extracted_root, capture_output=True, text=True)
    if res_cfg.returncode != 0:
        print("[FAIL] CMake configure failed in extracted archive:")
        print(res_cfg.stderr)
        return False

    # Build
    print("[*] Building project in extracted directory...")
    cmd_build = ["cmake", "--build", "build"]
    res_build = subprocess.run(cmd_build, cwd=extracted_root, capture_output=True, text=True)
    if res_build.returncode != 0:
        print("[FAIL] CMake build failed in extracted archive:")
        print(res_build.stderr)
        return False

    # Run CTest
    print("[*] Running CTest in extracted directory...")
    cmd_test = ["ctest", "--test-dir", "build", "--output-on-failure"]
    res_test = subprocess.run(cmd_test, cwd=extracted_root, capture_output=True, text=True)
    if res_test.returncode != 0:
        print("[FAIL] CTest failed in extracted archive:")
        print(res_test.stdout)
        print(res_test.stderr)
        return False

    print("[SUCCESS] All CTest targets passed in extracted source package build!")
    print(res_test.stdout.strip())

    # Run CLI demo
    cmd_demo = [os.path.join(extracted_root, "build", "thermo_cli.exe"), "--demo", "1"]
    res_demo = subprocess.run(cmd_demo, cwd=extracted_root, capture_output=True, text=True)
    if res_demo.returncode != 0 or "0.901844" not in res_demo.stdout:
        print("[FAIL] Demo test failed in extracted archive:")
        return False
    print("[SUCCESS] thermo_cli --demo 1 passed in extracted source package!")

    # Cleanup temporary directory
    shutil.rmtree(temp_dir, ignore_errors=True)
    print("[*] Cleaned up temporary extraction directory.")
    return True

if __name__ == "__main__":
    archive = create_source_archive()
    success = verify_extracted_build(archive)
    if not success:
        sys.exit(1)
    print("\n[ALL CHECKS PASSED] Source distribution archive is 100% verified!")
