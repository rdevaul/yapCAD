/**
 * Vite plugin to serve local yapCAD packages.
 * Handles /api/local-package?path=/path/to/package.ycpkg
 * Zips directories on the fly if needed.
 */
import { Plugin } from 'vite';
import { createReadStream, statSync, readdirSync, readFileSync } from 'fs';
import { join, resolve } from 'path';
import archiver from 'archiver';
import { PassThrough } from 'stream';

export function localPackagesPlugin(): Plugin {
  return {
    name: 'local-packages',
    configureServer(server) {
      server.middlewares.use('/api/local-package', async (req, res) => {
        const url = new URL(req.url || '', 'http://localhost');
        let pkgPath = url.searchParams.get('path');
        
        if (!pkgPath) {
          res.statusCode = 400;
          res.end(JSON.stringify({ error: 'Missing path parameter' }));
          return;
        }
        
        // Expand ~ to home directory
        if (pkgPath.startsWith('~')) {
          pkgPath = pkgPath.replace('~', process.env.HOME || '/');
        }
        
        // Resolve to absolute path
        pkgPath = resolve(pkgPath);
        
        console.log('[local-packages] Serving:', pkgPath);
        
        try {
          const stat = statSync(pkgPath);
          
          if (stat.isFile()) {
            // It's already a ZIP file, serve directly
            res.setHeader('Content-Type', 'application/zip');
            res.setHeader('Content-Length', stat.size);
            createReadStream(pkgPath).pipe(res);
          } else if (stat.isDirectory()) {
            // It's a directory, zip it on the fly
            res.setHeader('Content-Type', 'application/zip');
            
            const archive = archiver('zip', { zlib: { level: 5 } });
            
            archive.on('error', (err) => {
              console.error('[local-packages] Archive error:', err);
              res.statusCode = 500;
              res.end(JSON.stringify({ error: err.message }));
            });
            
            archive.pipe(res);
            
            // Add all files in the directory
            archive.directory(pkgPath, false);
            
            await archive.finalize();
          } else {
            res.statusCode = 400;
            res.end(JSON.stringify({ error: 'Path is neither a file nor directory' }));
          }
        } catch (err) {
          console.error('[local-packages] Error:', err);
          res.statusCode = 404;
          res.end(JSON.stringify({ error: `File not found: ${pkgPath}` }));
        }
      });
    },
  };
}

export default localPackagesPlugin;
