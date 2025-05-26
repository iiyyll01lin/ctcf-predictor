#!/usr/bin/env python3

"""
Proxy Detection and Environment Configuration for CTCF Predictor
This module provides dynamic proxy detection and configuration for Python scripts.
"""

import socket
import os
import subprocess
import json
from contextlib import closing
from typing import Dict, Optional, Tuple

class ProxyDetector:
    """Detects proxy availability and configures environment accordingly."""
    
    def __init__(self, proxy_host: str = "10.6.254.210", proxy_port: int = 3128, timeout: int = 5):
        """
        Initialize proxy detector.
        
        Args:
            proxy_host: Proxy server hostname or IP
            proxy_port: Proxy server port
            timeout: Connection timeout in seconds
        """
        self.proxy_host = proxy_host
        self.proxy_port = proxy_port
        self.timeout = timeout
        self.proxy_url = f"http://{proxy_host}:{proxy_port}"
    
    def check_proxy_connectivity(self) -> bool:
        """
        Check if proxy server is accessible.
        
        Returns:
            True if proxy is accessible, False otherwise
        """
        try:
            with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
                sock.settimeout(self.timeout)
                result = sock.connect_ex((self.proxy_host, self.proxy_port))
                return result == 0
        except Exception:
            return False
    
    def configure_environment(self) -> bool:
        """
        Configure environment based on proxy availability.
        
        Returns:
            True if proxy is available and configured, False if using direct connection
        """
        if self.check_proxy_connectivity():
            print(f"✓ Proxy accessible - configuring proxy settings ({self.proxy_url})")
            
            # Set proxy environment variables
            proxy_vars = {
                'HTTP_PROXY': self.proxy_url,
                'HTTPS_PROXY': self.proxy_url,
                'http_proxy': self.proxy_url,
                'https_proxy': self.proxy_url,
                'NO_PROXY': 'localhost,127.0.0.1,::1',
                'no_proxy': 'localhost,127.0.0.1,::1'
            }
            
            for var, value in proxy_vars.items():
                os.environ[var] = value
            
            return True
        else:
            print("✗ Proxy not accessible - using direct connection")
            # Clear proxy environment variables
            proxy_vars = ['HTTP_PROXY', 'HTTPS_PROXY', 'http_proxy', 'https_proxy']
            for var in proxy_vars:
                os.environ.pop(var, None)
            
            return False
    
    def get_requests_proxies(self) -> Optional[Dict[str, str]]:
        """
        Get proxy configuration for requests library.
        
        Returns:
            Proxy configuration dict for requests, or None for direct connection
        """
        if self.check_proxy_connectivity():
            return {
                'http': self.proxy_url,
                'https': self.proxy_url
            }
        return None
    
    def configure_pip_proxy(self) -> Tuple[bool, str]:
        """
        Configure pip to use proxy if available.
        
        Returns:
            Tuple of (proxy_available, pip_command_prefix)
        """
        if self.check_proxy_connectivity():
            return True, f"--proxy {self.proxy_url}"
        return False, ""
    
    def get_status(self) -> Dict[str, any]:
        """
        Get current proxy status and configuration.
        
        Returns:
            Dictionary with proxy status information
        """
        proxy_available = self.check_proxy_connectivity()
        
        return {
            "proxy_available": proxy_available,
            "proxy_host": self.proxy_host,
            "proxy_port": self.proxy_port,
            "proxy_url": self.proxy_url if proxy_available else None,
            "environment_configured": proxy_available,
            "connection_type": "proxy" if proxy_available else "direct"
        }

def main():
    """Main function for command-line usage."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Proxy detection and configuration")
    parser.add_argument("--host", default="10.6.254.210", help="Proxy host")
    parser.add_argument("--port", type=int, default=3128, help="Proxy port")
    parser.add_argument("--timeout", type=int, default=5, help="Connection timeout")
    parser.add_argument("--status", action="store_true", help="Show proxy status")
    parser.add_argument("--configure", action="store_true", help="Configure environment")
    
    args = parser.parse_args()
    
    detector = ProxyDetector(args.host, args.port, args.timeout)
    
    if args.status:
        status = detector.get_status()
        print(json.dumps(status, indent=2))
    elif args.configure:
        proxy_available = detector.configure_environment()
        exit(0 if proxy_available else 1)
    else:
        # Default behavior: check connectivity
        if detector.check_proxy_connectivity():
            print(f"✓ Proxy {args.host}:{args.port} is accessible")
            exit(0)
        else:
            print(f"✗ Proxy {args.host}:{args.port} is not accessible")
            exit(1)

if __name__ == "__main__":
    main()
