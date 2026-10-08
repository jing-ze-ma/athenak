# NOTE for DeltaAI: He giant bundle arrives as ONE TAR (start TASK-2026-10-07-deltaai-hegiant.md)

The user is uploading the bundle to DeltaAI now as a single tar instead of the directory copy described in section 3
of `TASK-2026-10-07-deltaai-hegiant.md`. Everything else in that TASK is unchanged: do it now.

- File: `hegiant_deltaai_bundle.tar` (2,277,775,360 bytes), md5 **6acf7ec374bb90fb98deb18f46dbe6e1**.
  It unpacks to `hegiant_deltaai_bundle/` (same 21 files + MD5SUMS as on viper).
- Ask the user where they put it if it is not in `/work/nvme/bivj/jma20/`. Then, in that directory:

```
echo "6acf7ec374bb90fb98deb18f46dbe6e1  hegiant_deltaai_bundle.tar" | md5sum -c
tar xf hegiant_deltaai_bundle.tar
bash $PWD/hegiant_deltaai_bundle/SETUP.sh $PWD/hegiant_deltaai_bundle
```

- If the md5 fails, the upload is incomplete: tell the user (resume with `rsync --partial`), do not unpack.
- After SETUP.sh passes, keep the bundle directory where it is for the whole run (TASK section 3); the tar may be
  deleted.
- Then: build f3a66907 with PROBLEM=he_star_m1 (TASK section 2), the 60-cycle smoke against viper job 12120235
  (TASK section 5), and start the 1-node x 4-GH200 chain only if the smoke passes.
- The viper chain hegiant897 (12119930-39) is still PENDING (no start estimate); leave it alone, the user decides
  which copy continues. Report the smoke numbers and the chain job ids in a NOTE on this branch.
