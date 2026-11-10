use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

pub const MAX_RESERVE_ATTEMPTS: u32 = 1000;

/// Creates the first unused candidate directory, giving up after a fixed
/// number of collisions so a misbehaving filesystem cannot hang tests.
pub fn reserve_directory(mut candidate: impl FnMut() -> PathBuf) -> std::io::Result<PathBuf> {
    for _ in 0..MAX_RESERVE_ATTEMPTS {
        let path = candidate();
        match fs::create_dir(&path) {
            Ok(()) => return Ok(path),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(e) => return Err(e),
        }
    }
    Err(std::io::Error::new(
        std::io::ErrorKind::AlreadyExists,
        format!("no unused directory name after {MAX_RESERVE_ATTEMPTS} attempts"),
    ))
}

/// Isolates all artifacts, including retained backups and directories, under one
/// collision-resistant root. Cleanup is best effort, including during unwinding.
pub struct TestDirectory(PathBuf);

impl TestDirectory {
    pub fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = reserve_directory(|| {
            std::env::temp_dir().join(format!(
                "aiips-test-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ))
        })
        .expect("create isolated test workspace");
        Self(path)
    }
}

impl AsRef<Path> for TestDirectory {
    fn as_ref(&self) -> &Path {
        &self.0
    }
}

impl std::ops::Deref for TestDirectory {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.0
    }
}

impl Drop for TestDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
