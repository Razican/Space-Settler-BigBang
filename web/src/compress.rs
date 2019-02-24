//! Response compression module.

use std::{
    io::{self, Cursor, Read},
    path::Path,
};

use failure::Error;
use flate2::{read::GzEncoder, Compression};
use maud::Markup;
use rocket::{
    http::{ContentType, Status},
    response::{NamedFile, Responder},
    Request, Response,
};

/// Compressed template.
#[derive(Debug)]
pub struct CompressedTemplate {
    /// Template to compress.
    template: Markup,
}

impl CompressedTemplate {
    /// Creates a new compressed template from a template.
    pub fn new(template: Markup) -> Self {
        Self { template }
    }
}

impl<'r> Responder<'r> for CompressedTemplate {
    fn respond_to(self, request: &Request) -> Result<Response<'r>, Status> {
        let mut response = self.template.respond_to(request)?;

        let headers = request.headers();
        // Check if requests accepts gzip encoding.
        if headers.contains("Accept")
            && headers.get("Accept-Encoding").any(
                |e| e.to_lowercase() == "gzip",
                // We compress the response here.
            )
            && compress_response(&mut response).is_err()
        {
            // Return an internal server error if compression went wrong.
            return Err(Status::InternalServerError);
        }
        Ok(response)
    }
}

/// Compressed file. Do not use for images.
#[derive(Debug)]
pub struct CompressedFile {
    /// File to compress.
    file: NamedFile,
    /// Content type for the response.
    content_type: ContentType,
}

impl CompressedFile {
    /// Creates a new compressed file from a path and a content type.
    pub fn new<P: AsRef<Path>>(file: P, content_type: ContentType) -> Result<Self, Error> {
        Ok(Self {
            file: NamedFile::open(file)?,
            content_type,
        })
    }
}

impl<'r> Responder<'r> for CompressedFile {
    fn respond_to(self, request: &Request) -> Result<Response<'r>, Status> {
        let mut response = self.file.respond_to(request)?;
        let _ = response.set_header(self.content_type);

        let headers = request.headers();
        // Check if requests accepts gzip encoding.
        if headers.contains("Accept")
            && headers.get("Accept-Encoding").any(
                |e| e.to_lowercase() == "gzip",
                // We compress the response here.
            )
            && compress_response(&mut response).is_err()
        {
            // Return an internal server error if compression went wrong.
            return Err(Status::InternalServerError);
        }
        Ok(response)
    }
}

/// Compresses the given response using Gzip.
///
/// Note that you should check if the client accepts compressed responses before compressing it.
fn compress_response(response: &mut Response) -> Result<(), io::Error> {
    let current_body = response.body_bytes().unwrap_or_default();
    let mut encoder = GzEncoder::new(&current_body[..], Compression::fast());

    let mut ret_body = Vec::new();
    let _ = encoder.read_to_end(&mut ret_body)?;

    let _ = response.set_raw_header("Content-Encoding", "gzip");
    response.set_sized_body(Cursor::new(ret_body));

    Ok(())
}
